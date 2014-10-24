/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once
#include "./qmc_data.hpp"

namespace cthyb {

// Move a C or C^dagger operator to a different time
class move_shift_operator {

 qmc_data& data;
 configuration& config;
 mc_tools::random_generator& rng;
 bool record_histograms;
 std::map<std::string, statistics::histogram_segment_bin> histos; // Analysis histograms
 double delta_tau;
 qmc_data::trace_t new_trace;
 time_pt tau_old, tau_new;
 op_desc op_old, op_new;
 typedef det_manip::det_manip<qmc_data::delta_block_adaptor> det_type;
 det_type::RollDirection roll_direction;
 int block_index;

 public:
 //-----------------------------------------------

 move_shift_operator(qmc_data& data, mc_tools::random_generator& rng, bool record_histograms)
    : data(data),
      config(data.config),
      rng(rng),
      record_histograms(record_histograms) {
  if (record_histograms) { 
   //FIXME grid too coarse, need to update
   histos.insert({"length_proposed", {0, config.beta(), 100, "hist_length_proposed.dat"}});
   histos.insert({"length_accepted", {0, config.beta(), 100, "hist_length_accepted.dat"}});
  }
 }

 //---------------------

 mc_weight_type attempt() {

// FIXME to be studied: effect of creating move_shift per block vs "global" move_shift
//  // --- Choose an operator in block to move at random and whether dagger or not
//  auto& det = data.dets[block_index];
//  int det_size = det.size();
//  if (det_size == 0) return 0; // nothing to remove
//  int num_op_old = rng(det_size);
//  bool is_dagger = rng(0,1);

  // --- Choose an operator in configuration to shift at random
  // By choosing an *operator* in config directly, not bias based on det size introduced
  const int config_size(config.size());
  if (config_size==0) return 0;
  const int op_pos_in_config((rng(config_size)));

  // --- Find operator (and its characteristics) from the configuration
  auto itconfig = config.begin();
  for (int i=0; i<op_pos_in_config; ++i, ++itconfig) {} // Get to right position in config
  tau_old = (*itconfig).first;
  op_old = (*itconfig).second;
  block_index = op_old.block_index;
  auto is_dagger = op_old.dagger;
  // Properties corresposing to det
  auto& det = data.dets[block_index];
  int det_size = det.size();
  if (det_size == 0) return 0; // nothing to remove

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "* Attempt for move_shift_operator (block " << block_index << ")" << std::endl;
  std::cerr << "* Configuration before:" << std::endl;
  std::cerr << config;
  data.imp_trace.tree.graphviz(std::ofstream("tree_before"));
#endif

  // Construct new operator
  // Choose a new inner index (this is done here for compatibility) -- FIXME alpha needs to change?
  auto inner_new = rng(data.n_inner[block_index]);
  op_new = op_desc{block_index, inner_new, is_dagger, data.linindex[std::make_pair(block_index, inner_new)]};

  // --- Determine new time to shift the operator to.
  // The time must fall in the range between the closest operators on the left and 
  // right belonging to the same block. First determine these.

  time_pt tR, tL;
  int ic_dag = 0, ic = 0;
  int op_pos_in_det;

  if (det_size > 1) {

    // Get the position op_pos_in_det of op_old in the det.
    // Find the c and c_dag operators at the right of op_old (at smaller times)
    // They could be the last entries (earliest times)

    for (ic_dag = 0; ic_dag < det_size; ++ic_dag) { // c_dag
     if (det.get_x(ic_dag).first < tau_old) break;
    } 
    for (ic = 0; ic < det_size; ++ic) { // c
     if (det.get_y(ic).first < tau_old) break;
    } 

    op_pos_in_det = (is_dagger ? ic_dag : ic); // This finds the operator on the right
    --op_pos_in_det; // Rewind by one to find the operator

    // Find the times of the operator at the right of op_old with cyclicity
    auto tRdag   = ( ic_dag != det_size ? det.get_x(ic_dag).first : det.get_x(0).first );
    auto tRnodag =     ( ic != det_size ? det.get_y(ic).first     : det.get_y(0).first );
    // Then deduce the closest one and put its distance to op_old in tR
    tR = ((tau_old - tRdag) > (tau_old - tRnodag) ? tRnodag : tRdag);

    // Reset iterator to op_old position
    if (is_dagger) --ic_dag; else --ic; 

    // Find the times of the operator at the left of op_old with cyclicity
    auto tLdag   = ( ic_dag != 0 ? det.get_x(--ic_dag).first : det.get_x(det_size-1).first );
    auto tLnodag =     ( ic != 0 ? det.get_y(--ic).first     : det.get_y(det_size-1).first );
    // Then deduce the closest one and put its distance to op_old in tL
    tL = ((tLdag - tau_old) > (tLnodag - tau_old) ? tLnodag : tLdag);
    // Choose new random time
    tau_new = tR + time_pt::random(rng, tL - tR, config.beta());

  } else { // det_size = 1

    op_pos_in_det = 0;
    // Choose new random time, can be anywhere between beta and 0
    tau_new = time_pt::random(rng, config.beta(), config.beta());

  }

  // Record the length of the proposed shift
  delta_tau = double(tau_new - tau_old);
  if (record_histograms) histos["length_proposed"] << delta_tau;

  // --- Modify the tree
  // Mark the operator at original time for deletion in the tree
  auto tau_old_tree = data.imp_trace.try_delete(op_pos_in_det, block_index, is_dagger);
  // FIXME could check if tau_old_tree == tau_old

#ifdef EXT_DEBUG
  std::cerr << "* Proposing to shift:" << std::endl;
  std::cerr << (is_dagger ? "Cdag" : "C");
  std::cerr << "(" << op_old.block_index << "," << op_old.inner_index << ")";
  std::cerr << " tau = " << tau_old << std::endl;
  std::cerr << " to " << std::endl;
  std::cerr << (is_dagger ? "Cdag" : "C");
  std::cerr << "(" << op_new.block_index << "," << op_new.inner_index << ")";
  std::cerr << " tau = " << tau_new << std::endl;
#endif

  // Try to insert the new operator at shifted time in the tree
  try {
   data.imp_trace.try_insert(tau_new, op_new);
  }
  catch (rbt_insert_error const&) {
   //std::cerr << "Insert error : recovering ... " << std::endl;
   data.imp_trace.cancel_insert();
   return 0;
  }

  // --- Compute the det ratio

  // Do we need to roll the determinant?
  typedef det_manip::det_manip<qmc_data::delta_block_adaptor> det_type;
  roll_direction = det_type::None;

  // Check if we went through \tau = 0 or \tau = \beta
  if (det_size > 1) {
   if ((tau_old > tL) && (tau_new < tL)) roll_direction = (is_dagger ? det_type::Up   : det_type::Left);
   if ((tau_old < tL) && (tau_new > tL)) roll_direction = (is_dagger ? det_type::Down : det_type::Right); 
  }

  // Replace old row/column with new operator time/inner_index. Returns the ratio of dets (Cf det_manip doc).
  auto det_ratio = (is_dagger ? det.try_change_row(op_pos_in_det, {tau_new,op_new.inner_index}) : det.try_change_col(op_pos_in_det, {tau_new,op_new.inner_index}));

  // for quick abandon
  double random_number = rng.preview();
  if (random_number == 0.0) return 0;
  double p_yee = std::abs(det_ratio / data.trace);

  // --- Compute the trace ratio
  new_trace = data.imp_trace.estimate(p_yee, random_number);
  if (new_trace == 0.0) {
#ifdef EXT_DEBUG
   std::cout << "trace == 0" << std::endl;
#endif
   return 0;
  }
  auto trace_ratio = new_trace / data.trace;
  if (!std::isfinite(trace_ratio)) TRIQS_RUNTIME_ERROR << "trace_ratio not finite" << new_trace << "  "<< data.trace<<"  "<< new_trace /data.trace ;

  // --- Compute the weight
  mc_weight_type p = trace_ratio * det_ratio;

#ifdef EXT_DEBUG
  std::cerr << "Trace ratio: " << trace_ratio << '\t';
  std::cerr << "Det ratio: " << det_ratio << '\t';
  std::cerr << "Weight: " << p << std::endl;
  std::cerr << "p_yee* newtrace: " << p_yee * new_trace<< std::endl;
#endif

  return p;
 }

 //----------------

 mc_weight_type accept() {

  // Update the tree
  data.imp_trace.confirm_shift();

  // Update the configuration 
  config.erase(tau_old);
  config.insert(tau_new, op_new);
 
  // Update the determinant
  data.dets[block_index].complete_operation();
  data.update_sign();
  data.trace = new_trace;
  if (record_histograms) histos["length_accepted"] << delta_tau;

#ifdef EXT_DEBUG
  std::cerr << "* Configuration after: " << std::endl;
  std::cerr << config;
#endif

  return data.current_sign / data.old_sign * data.dets[block_index].roll_matrix(roll_direction);
 }

 //----------------

 void reject() {
  data.imp_trace.cancel_shift();
#ifdef EXT_DEBUG
  std::cerr << "* Configuration after: " << std::endl;
  std::cerr << config;
#endif
 }
};
}
