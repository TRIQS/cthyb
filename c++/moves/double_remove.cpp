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

#include "./double_remove.hpp"

namespace cthyb {

// -----------------------------------------------------------------------------
histogram* move_remove_c_c_cdag_cdag::add_histo(std::string const& name, histo_map_t* histos) {
 if (!histos) return nullptr;
 auto new_histo = histos->insert({name, {.0, config.beta(), 100}});
 return &(new_histo.first->second);
}

// -----------------------------------------------------------------------------
mc_weight_t move_remove_c_c_cdag_cdag::attempt() {

#ifdef EXT_DEBUG
 std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
 std::cerr << "In config " << config.get_id() << std::endl;
 std::cerr << "* Attempt for move_remove_c_c_cdag_cdag (blocks " << block_index1 << ", " << block_index2 << ")" << std::endl;
#endif

 auto& det1 = data.dets[block_index1];
 auto& det2 = data.dets[block_index2];
 det_scalar_t det_ratio;

 // Pick two pairs of C, Cdagger to remove at random
 // Remove the operators from the traces
 int det1_size = det1.size();
 int det2_size = det2.size();
 if (block_index1 == block_index2) {
  if (det1_size < 2) return 0; // det does not contain two operators
 } else {
  if ((det1_size == 0) || (det2_size == 0)) return 0; // one of two dets is empty
 }
 int num_c_dag1 = rng(det1_size), num_c1 = rng(det1_size);
 int num_c_dag2 = rng(det2_size), num_c2 = rng(det2_size);
 if ((block_index1 == block_index2) && ((num_c_dag1 == num_c_dag2) || (num_c1 == num_c2)))
  return 0; // picked the same operator twice

#ifdef EXT_DEBUG
 std::cerr << "* Proposing to remove: ";
 std::cerr << num_c_dag1 << "-th Cdag(" << block_index1 << ",...), ";
 std::cerr << num_c1 << "-th C(" << block_index1 << ",...)" << std::endl;
 std::cerr << " and ";
 std::cerr << num_c_dag2 << "-th Cdag(" << block_index2 << ",...), ";
 std::cerr << num_c2 << "-th C(" << block_index2 << ",...)" << std::endl;
#endif

 // now mark 2 nodes for deletion
 tau1 = data.imp_trace.try_delete(num_c1, block_index1, false);
 tau2 = data.imp_trace.try_delete(num_c_dag1, block_index1, true);
 tau3 = data.imp_trace.try_delete(num_c2, block_index2, false);
 tau4 = data.imp_trace.try_delete(num_c_dag2, block_index2, true);

 dtau1 = double(tau2 - tau1);
 dtau2 = double(tau4 - tau3);
 if (histo_proposed1) {
  *histo_proposed1 << dtau1;
  *histo_proposed2 << dtau2;
 }

 if (block_index1 == block_index2) {
  det_ratio = det1.try_remove2(num_c_dag1, num_c_dag2, num_c1, num_c2);
 } else { // block_index1 != block_index2
  auto det_ratio1 = det1.try_remove(num_c_dag1, num_c1);
  auto det_ratio2 = det2.try_remove(num_c_dag2, num_c2);
  det_ratio = det_ratio1 * det_ratio2;
 }

 // proposition probability
 mc_weight_t t_ratio;
 // Note: Must use the size of the det before the try_delete!
 if (block_index1 == block_index2) {
  // Here, we use the fact that the two cdag/c proposed to be removed in the det can be at the same
  // positions in the det, and thus remove prob is NOT (detsize+2)*(detsize+1)
  t_ratio = std::pow(block_size1 * config.beta() / double(det1_size), 4);
 } else {
  t_ratio =
      std::pow(block_size1 * config.beta() / double(det1_size), 2) * std::pow(block_size2 * config.beta() / double(det2_size), 2);
 }

 // For quick abandon
 double random_number = rng.preview();
 if (random_number == 0.0) return 0;
 double p_yee = std::abs(det_ratio / t_ratio / data.atomic_weight);

 // recompute the trace
 std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
 if (new_atomic_weight == 0.0) {
#ifdef EXT_DEBUG
  std::cerr << "atomic_weight == 0" << std::endl;
#endif
  return 0;
 }
 auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
 if (!isfinite(atomic_weight_ratio))
  TRIQS_RUNTIME_ERROR << "atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                      << new_atomic_weight / data.atomic_weight << " in config " << config.get_id();

 mc_weight_t p = atomic_weight_ratio * det_ratio;

#ifdef EXT_DEBUG
 std::cerr << "Trace ratio: " << atomic_weight_ratio << '\t';
 std::cerr << "Det ratio: " << det_ratio << '\t';
 std::cerr << "Prefactor: " << t_ratio << '\t';
 std::cerr << "Weight: " << p / t_ratio << std::endl;
#endif

 if (!isfinite(p)) TRIQS_RUNTIME_ERROR << "(remove) p not finite :" << p << " in config " << config.get_id();
 if (!isfinite(p / t_ratio))
  TRIQS_RUNTIME_ERROR << "p / t_ratio not finite p : " << p << " t_ratio :  " << t_ratio << " in config " << config.get_id();
 return p / t_ratio;
}

// -----------------------------------------------------------------------------
mc_weight_t move_remove_c_c_cdag_cdag::accept() {

 // remove from the tree
 data.imp_trace.confirm_delete();

 // remove from the configuration
 config.erase(tau1);
 config.erase(tau2);
 config.erase(tau3);
 config.erase(tau4);
 config.finalize();

 // remove from the determinants
 if (block_index1 == block_index2) {
  data.dets[block_index1].complete_operation();
 } else {
  data.dets[block_index1].complete_operation();
  data.dets[block_index2].complete_operation();
 }
 data.update_sign();

 data.atomic_weight = new_atomic_weight;
 data.atomic_reweighting = new_atomic_reweighting;

 if (histo_accepted1) {
  *histo_accepted1 << dtau1;
  *histo_accepted2 << dtau2;
 }

#ifdef EXT_DEBUG
 std::cerr << "* Move move_remove_c_c_cdag_cdag accepted" << std::endl;
 std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
 check_det_sequence(data.dets[block_index1], config.get_id());
 check_det_sequence(data.dets[block_index2], config.get_id());
#endif

 return data.current_sign / data.old_sign;
}

// -----------------------------------------------------------------------------
void move_remove_c_c_cdag_cdag::reject() {

 config.finalize();
 data.imp_trace.cancel_delete();

#ifdef EXT_DEBUG
 std::cerr << "* Move move_remove_c_c_cdag_cdag rejected" << std::endl;
 std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
 check_det_sequence(data.dets[block_index1], config.get_id());
 check_det_sequence(data.dets[block_index2], config.get_id());
#endif
}
}
