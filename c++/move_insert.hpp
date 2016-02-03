/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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

// Insertion of C, C^dagger operator
class move_insert_c_cdag {

 qmc_data& data;
 configuration& config;
 mc_tools::random_generator& rng;
 int block_index, block_size;
 bool performance_analysis;
 std::map<std::string, statistics::histogram_segment_bin> histos; // Analysis histograms
 double dtau;
 h_scalar_t new_atomic_weight, new_atomic_reweighting;
 time_pt tau1, tau2;
 op_desc op1, op2;

 public:
 //-----------------------------------------------

 move_insert_c_cdag(int block_index, int block_size, qmc_data& data, mc_tools::random_generator& rng, bool performance_analysis)
    : data(data),
      config(data.config),
      rng(rng),
      block_index(block_index),
      block_size(block_size),
      performance_analysis(performance_analysis) {
  if (performance_analysis) {
   histos.insert({"insert_length_proposed", {0, config.beta(), 100, "histo_insert_length_proposed.dat"}});
   histos.insert({"insert_length_accepted", {0, config.beta(), 100, "histo_insert_length_accepted.dat"}});
  }
 }

 //---------------------

 mc_weight_t attempt() {

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "In config " << config.get_id() << std::endl;
  std::cerr << "* Attempt for move_insert_c_cdag (block " << block_index << ")" << std::endl;
#endif

  // Pick up the value of alpha and choose the operators
  auto rs1 = rng(block_size), rs2 = rng(block_size);
  op1 = op_desc{block_index, rs1, true, data.linindex[std::make_pair(block_index, rs1)]};
  op2 = op_desc{block_index, rs2, false, data.linindex[std::make_pair(block_index, rs2)]};

  // Choice of times for insertion. Find the time as double and them put them on the grid.
  tau1 = data.tau_seg.get_random_pt(rng);
  tau2 = data.tau_seg.get_random_pt(rng);

#ifdef EXT_DEBUG
  std::cerr << "* Proposing to insert:" << std::endl;
  std::cerr << op1 << " at " << tau1 << std::endl;
  std::cerr << op2 << " at " << tau2 << std::endl;
#endif

  // record the length of the proposed insertion
  dtau = double(tau2 - tau1);
  if (performance_analysis) histos["insert_length_proposed"] << dtau;

  // Insert the operators op1 and op2 at time tau1, tau2
  // 1- In the very exceptional case where the insert has failed because an operator is already sitting here
  // (cf std::map doc for insert return), we reject the move.
  // 2- If ok, we store the iterator to the inserted operators for later removal in reject if necessary
  try {
   data.imp_trace.try_insert(tau1, op1);
   data.imp_trace.try_insert(tau2, op2);
  }
  catch (rbt_insert_error const&) {
   std::cerr << "Insert error : recovering ... " << std::endl;
   data.imp_trace.cancel_insert();
   return 0;
  }
 
  // Computation of det ratio
  auto& det = data.dets[block_index];
  int det_size = det.size();

  // Find the position for insertion in the determinant
  // NB : the determinant stores the C in decreasing time order.
  int num_c_dag, num_c;
  for (num_c_dag = 0; num_c_dag < det_size; ++num_c_dag) {
   if (det.get_x(num_c_dag).first < tau1) break;
  }
  for (num_c = 0; num_c < det_size; ++num_c) {
   if (det.get_y(num_c).first < tau2) break;
  }

  // Insert in the det. Returns the ratio of dets (Cf det_manip doc).
  auto det_ratio = det.try_insert(num_c_dag, num_c, {tau1, op1.inner_index}, {tau2, op2.inner_index});

  // proposition probability
  mc_weight_t t_ratio = std::pow(block_size * config.beta() / double(det.size() + 1), 2);

  // For quick abandon
  double random_number = rng.preview();
  if (random_number == 0.0) return 0;
  double p_yee = std::abs(t_ratio * det_ratio / data.atomic_weight);

  // computation of the new trace after insertion
  std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
  if (new_atomic_weight == 0.0) {
#ifdef EXT_DEBUG
   std::cerr << "atomic_weight == 0" << std::endl;
#endif
   return 0;
  }
  auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
  if (!isfinite(atomic_weight_ratio))
   TRIQS_RUNTIME_ERROR << "trace_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                       << new_atomic_weight / data.atomic_weight << " in config " << config.get_id();

  mc_weight_t p = atomic_weight_ratio * det_ratio;

#ifdef EXT_DEBUG
  std::cerr << "Atomic ratio: " << atomic_weight_ratio << '\t';
  std::cerr << "Det ratio: " << det_ratio << '\t';
  std::cerr << "Prefactor: " << t_ratio << '\t';
  std::cerr << "Weight: " << p* t_ratio << std::endl;
  std::cerr << "p_yee * newtrace: " << p_yee * new_atomic_weight<< std::endl;
#endif

  if (!isfinite(p * t_ratio)) TRIQS_RUNTIME_ERROR << "p * t_ratio not finite p : " << p << " t_ratio : "<< t_ratio << " in config " << config.get_id();
  return p * t_ratio;
 }

 //----------------

 mc_weight_t accept() {

  // insert in the tree
  data.imp_trace.confirm_insert();

  // insert in the configuration
  config.insert(tau1, op1);
  config.insert(tau2, op2);
  config.finalize();

  // insert in the determinant
  data.dets[block_index].complete_operation();
  data.update_sign();
  data.atomic_weight = new_atomic_weight;
  data.atomic_reweighting = new_atomic_reweighting;
  if (performance_analysis) histos["insert_length_accepted"] << dtau;

#ifdef EXT_DEBUG
  std::cerr << "* Move move_insert_c_cdag accepted" << std::endl;
  std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  check_det_sequence(data.dets[block_index],config.get_id());
#endif

  return data.current_sign / data.old_sign;
 }

 //----------------

 void reject() {

  config.finalize();
  data.imp_trace.cancel_insert();

#ifdef EXT_DEBUG
  std::cerr << "* Move move_insert_c_cdag rejected" << std::endl;
  std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  check_det_sequence(data.dets[block_index],config.get_id());
#endif
 }
};
}
