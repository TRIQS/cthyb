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
class move_insert_c_c_cdag_cdag {

 qmc_data& data;
 configuration& config;
 mc_tools::random_generator& rng;
 int block_index1, block_index2, block_size1, block_size2;
 bool performance_analysis;
 std::map<std::string, statistics::histogram_segment_bin> histos; // Analysis histograms
 double delta_tau1, delta_tau2;
 double new_atomic_weight;
 qmc_data::trace_t new_atomic_reweighting;
 time_pt tau1, tau2, tau3, tau4;
 op_desc op1, op2, op3, op4;

 public:
 //-----------------------------------------------

 move_insert_c_c_cdag_cdag(int block_index1, int block_index2, int block_size1, int block_size2, qmc_data& data, mc_tools::random_generator& rng, bool performance_analysis)
    : data(data),
      config(data.config),
      rng(rng),
      block_index1(block_index1),
      block_size1(block_size1),
      block_index2(block_index2),
      block_size2(block_size2),
      performance_analysis(performance_analysis) {
  if (performance_analysis) {
   histos.insert({"double_insert_length_proposed", {0, config.beta(), 100, "histo_double_insert_length_proposed.dat"}});
   histos.insert({"double_insert_length_accepted", {0, config.beta(), 100, "histo_double_insert_length_accepted.dat"}});
  }
 }

 //---------------------

 mc_weight_type attempt() {

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "In config " << config.get_id() << std::endl;
  std::cerr << "* Attempt for move_insert_c_c_cdag_cdag (blocks " << block_index1 << ", " << block_index2 << ")" << std::endl;
#endif

  // Pick up the value of alpha and choose the operators
  auto rs1 = rng(block_size1), rs2 = rng(block_size1), rs3 = rng(block_size2), rs4 = rng(block_size2);
  op1 = op_desc{block_index1, rs1, true, data.linindex[std::make_pair(block_index1, rs1)]};
  op2 = op_desc{block_index1, rs2, false, data.linindex[std::make_pair(block_index1, rs2)]};
  op3 = op_desc{block_index2, rs3, true, data.linindex[std::make_pair(block_index2, rs3)]};
  op4 = op_desc{block_index2, rs4, false, data.linindex[std::make_pair(block_index2, rs4)]};

  // Choice of times for insertion. Find the time as double and them put them on the grid.
  tau1 = data.tau_seg.get_random_pt(rng);
  tau2 = data.tau_seg.get_random_pt(rng);
  tau3 = data.tau_seg.get_random_pt(rng);
  tau4 = data.tau_seg.get_random_pt(rng);
  if ((tau1 == tau3) or (tau2 == tau4)) return 0; // trying to insert/remove two operators at exactly the same time

#ifdef EXT_DEBUG
  std::cerr << "* Proposing to insert:" << std::endl;
  std::cerr << op1 << " at " << tau1 << std::endl;
  std::cerr << op2 << " at " << tau2 << std::endl;
  std::cerr << op3 << " at " << tau3 << std::endl;
  std::cerr << op4 << " at " << tau4 << std::endl;
#endif

  // record the length of the proposed insertion
  delta_tau1 = double(tau2 - tau1);
  delta_tau2 = double(tau4 - tau3);
  if (performance_analysis) {
   histos["double_insert_length_proposed"] << delta_tau1;
   histos["double_insert_length_proposed"] << delta_tau2;
  }

  // Insert the operators op1, op2, op3, op4 at time tau1, tau2, tau3, tau4
  // 1- In the very exceptional case where the insert has failed because an operator is already sitting here
  // (cf std::map doc for insert return), we reject the move.
  // 2- If ok, we store the iterator to the inserted operators for later removal in reject if necessary
  try {
   data.imp_trace.try_insert(tau1, op1);
   data.imp_trace.try_insert(tau2, op2);
   data.imp_trace.try_insert(tau3, op3);
   data.imp_trace.try_insert(tau4, op4);
  }
  catch (rbt_insert_error const&) {
   std::cerr << "Insert error : recovering ... " << std::endl;
   data.imp_trace.cancel_insert();
   return 0;
  }
 
  // Computation of det ratio
  auto& det1 = data.dets[block_index1];
  auto& det2 = data.dets[block_index2];
  int det1_size = det1.size();
  int det2_size = det2.size();
  double det_ratio;
  
  // Find the position for insertion in the determinant
  // NB : the determinant stores the C in decreasing time order.
  int num_c_dag1, num_c1, num_c_dag2, num_c2;
  for (num_c_dag1 = 0; num_c_dag1 < det1_size; ++num_c_dag1) {
   if (det1.get_x(num_c_dag1).first < tau1) break;
  }
  for (num_c1 = 0; num_c1 < det1_size; ++num_c1) {
   if (det1.get_y(num_c1).first < tau2) break;
  }
  for (num_c_dag2 = 0; num_c_dag2 < det2_size; ++num_c_dag2) {
   if (det2.get_x(num_c_dag2).first < tau3) break;
  }
  for (num_c2 = 0; num_c2 < det2_size; ++num_c2) {
   if (det2.get_y(num_c2).first < tau4) break;
  }
  
  // Insert in the det. Returns the ratio of dets (Cf det_manip doc).
  if (block_index1 == block_index2) {
   // The determinant positions that need to be passed to det_manip are those in the *final* det of size N+2.
   // Shift the operator at the smaller time one step further in the determinant to account for the larger operator. 
   // This shfit must be done in general, and not only when num_c(_dag)1 and num_c(dag_)2 are the same!!
   if (tau1 < tau3) num_c_dag1++; else num_c_dag2++;
   if (tau2 < tau4) num_c1++; else num_c2++;
   det_ratio = det1.try_insert2(num_c_dag1, num_c_dag2, num_c1, num_c2, {tau1, op1.inner_index}, {tau3, op3.inner_index}, 
                                                                             {tau2, op2.inner_index}, {tau4, op4.inner_index});
  } else {
   auto det_ratio1 = det1.try_insert(num_c_dag1, num_c1, {tau1, op1.inner_index}, {tau2, op2.inner_index});
   auto det_ratio2 = det2.try_insert(num_c_dag2, num_c2, {tau3, op3.inner_index}, {tau4, op4.inner_index});
   det_ratio = det_ratio1 * det_ratio2;
  }

  // proposition probability
  double t_ratio; 
  if (block_index1 == block_index2) {
   // (ways to insert 4 operators in det1)/((ways to remove 4 operators from det that is larger by two))
   // Here, we use the fact that the two cdag/c proposed to be removed in the det can be at the same 
   // positions in the det, and thus remove prob is NOT (detsize+2)*(detsize+1)
   t_ratio = std::pow(block_size1 * config.beta() / double(det1.size() + 2), 4);
  } else { 
   // product of two separate inserts, one in det1 and one in det2
   t_ratio = std::pow(block_size1 * config.beta() / double(det1.size() + 1), 2) * std::pow(block_size2 * config.beta() / double(det2.size() + 1), 2);
  }

  // For quick abandon
  double random_number = rng.preview();
  if (random_number == 0.0) return 0;
  double p_yee = std::abs(t_ratio * det_ratio / data.atomic_weight);

  // computation of the new atomic_weight after insertion
  std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
  if (new_atomic_weight == 0.0) {
#ifdef EXT_DEBUG
   std::cerr << "atomic_weight == 0" << std::endl;
#endif
   return 0;
  }
  auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
  if (!std::isfinite(atomic_weight_ratio)) TRIQS_RUNTIME_ERROR << "atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " " << new_atomic_weight /data.atomic_weight << " in config " << config.get_id();

  mc_weight_type p = atomic_weight_ratio * det_ratio;

#ifdef EXT_DEBUG
  std::cerr << "Trace ratio: " << atomic_weight_ratio << '\t';
  std::cerr << "Det ratio: " << det_ratio << '\t';
  std::cerr << "Prefactor: " << t_ratio << '\t';
  std::cerr << "Weight: " << p* t_ratio << std::endl;
  std::cerr << "p_yee * newtrace: " << p_yee * new_atomic_weight<< std::endl;
#endif

  if (!std::isfinite(p * t_ratio)) TRIQS_RUNTIME_ERROR << "p * t_ratio not finite p : " << p << " t_ratio :  "<< t_ratio << " in config " << config.get_id();
  return p * t_ratio;
 }

 //----------------

 mc_weight_type accept() {

  // insert in the tree
  data.imp_trace.confirm_insert();

  // insert in the configuration
  config.insert(tau1, op1);
  config.insert(tau2, op2);
  config.insert(tau3, op3);
  config.insert(tau4, op4);
  config.finalize();

  // insert in the determinant
  if (block_index1 == block_index2) {
   data.dets[block_index1].complete_operation();
  } else {
   data.dets[block_index1].complete_operation();
   data.dets[block_index2].complete_operation();
  }
  data.update_sign();

  data.atomic_weight = new_atomic_weight;
  data.atomic_reweighting = new_atomic_reweighting;

  if (performance_analysis) {
   histos["double_insert_length_accepted"] << delta_tau1;
   histos["double_insert_length_accepted"] << delta_tau2;
  }

#ifdef EXT_DEBUG
  std::cerr << "* Move move_insert_c_c_cdag_cdag accepted" << std::endl;
  std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  check_det_sequence(data.dets[block_index1],config.get_id());
  check_det_sequence(data.dets[block_index2],config.get_id());
#endif

  return data.current_sign / data.old_sign;
 }

 //----------------

 void reject() {

  config.finalize();
  data.imp_trace.cancel_insert();

#ifdef EXT_DEBUG
  std::cerr << "* Move move_insert_c_c_cdag_cdag rejected" << std::endl;
  std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
  check_det_sequence(data.dets[block_index1],config.get_id());
  check_det_sequence(data.dets[block_index2],config.get_id());
#endif

 }
};
}
