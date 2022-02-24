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
#include "./util.hpp"

namespace triqs_cthyb {

  move_remove_c_c_cdag_cdag::move_remove_c_c_cdag_cdag(int block_index1, int block_index2, int block_size1, int block_size2,
                                                       std::string const &block_name1, std::string const &block_name2, qmc_data &data,
                                                       mc_tools::random_generator &rng, histo_map_t *histos)
     : data(data),
       config(data.config),
       rng(rng),
       block_index1(block_index1),
       block_index2(block_index2),
       block_size1(block_size1),
       block_size2(block_size2),
       histo_proposed1(add_histo("double_remove_length_proposed_" + block_name1, histos, config.beta())),
       histo_proposed2(add_histo("double_remove_length_proposed_" + block_name2, histos, config.beta())),
       histo_accepted1(add_histo("double_remove_length_accepted_" + block_name1, histos, config.beta())),
       histo_accepted2(add_histo("double_remove_length_accepted_" + block_name2, histos, config.beta())) {}

  mc_weight_t move_remove_c_c_cdag_cdag::attempt() {

    LOG("{}\n* Attempt for move_remove_c_c_cdag_cdag (block {}, {})", debug_config_print_start, block_index1, block_index2);

    auto &det1 = data.dets[block_index1];
    auto &det2 = data.dets[block_index2];
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
    if ((block_index1 == block_index2) && ((num_c_dag1 == num_c_dag2) || (num_c1 == num_c2))) return 0; // picked the same operator twice

    LOG("* Proposing to remove:");
    LOG("    for block {}:\n       {}-th Cdag\n       {}-th C", block_index1, num_c_dag1, num_c1);
    LOG("    for block {}:\n       {}-th Cdag\n       {}-th C", block_index2, num_c_dag2, num_c2);

    // now mark 2 nodes for deletion
    tau1 = det1.get_y(num_c1).first;
    tau2 = det1.get_x(num_c_dag1).first;
    tau3 = det2.get_y(num_c2).first;
    tau4 = det2.get_x(num_c_dag2).first;

    data.imp_trace.try_delete(tau1);
    data.imp_trace.try_delete(tau2);
    data.imp_trace.try_delete(tau3);
    data.imp_trace.try_delete(tau4);

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
      det_ratio       = det_ratio1 * det_ratio2;
    }

    // proposition probability
    mc_weight_t t_ratio;
    // Note: Must use the size of the det before the try_delete!
    if (block_index1 == block_index2) {
      // Here, we use the fact that the two cdag/c proposed to be removed in the det can be at the same
      // positions in the det, and thus remove prob is NOT (detsize+2)*(detsize+1)
      t_ratio = std::pow(block_size1 * config.beta() / double(det1_size), 4);
    } else {
      t_ratio = std::pow(block_size1 * config.beta() / double(det1_size), 2) * std::pow(block_size2 * config.beta() / double(det2_size), 2);
    }

    // For quick abandon
    double random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double p_yee = std::abs(det_ratio / t_ratio / data.atomic_weight);

    // recompute the trace
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) {
      LOG("atomic_weight == 0");
      return 0;
    }
    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                          << new_atomic_weight / data.atomic_weight << " in config " << config.get_id();

    mc_weight_t p = atomic_weight_ratio * det_ratio;

    LOG("Trace ratio: {}  Det ratio: {}  Prefactor: {} ", atomic_weight_ratio, det_ratio, t_ratio);
    LOG("Weight: {}\n  p_yee * newtrace: {} ", p / t_ratio, p_yee * new_atomic_weight);
    ALWAYS_EXPECTS(isfinite(p), "(remove) p = {} not finite", p);
    ALWAYS_EXPECTS(isfinite(p), "(remove) p / t_ratio not finite. p = {}, t_ratio = {}, config :\n {} ", p, t_ratio, config);

    return p / t_ratio;
  }

  // -----------------------------------------------------------

  mc_weight_t move_remove_c_c_cdag_cdag::accept() {

    // remove from the tree
    data.imp_trace.confirm_delete();

    // remove from the configuration
    config.erase(tau1);
    config.erase(tau2);
    config.erase(tau3);
    config.erase(tau4);

    // remove from the determinants
    if (block_index1 == block_index2) {
      data.dets[block_index1].complete_operation();
    } else {
      data.dets[block_index1].complete_operation();
      data.dets[block_index2].complete_operation();
    }
    data.update_sign();

    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;

    if (histo_accepted1) {
      *histo_accepted1 << dtau1;
      *histo_accepted2 << dtau2;
    }

    LOG("* Accepted \n{}", debug_config_print_end);
#ifdef EXT_DEBUG
    check_det_sequence(data.dets[block_index1], config.get_id());
    check_det_sequence(data.dets[block_index2], config.get_id());
#endif

    return data.current_sign / data.old_sign;
  }

  // -----------------------------------------------------------

  void move_remove_c_c_cdag_cdag::reject() {

    data.imp_trace.cancel_delete();
    // remove from the determinants
    if (block_index1 == block_index2) {
      data.dets[block_index1].reject_last_try();
    } else {
      data.dets[block_index1].reject_last_try();
      data.dets[block_index2].reject_last_try();
    }

#ifdef EXT_DEBUG
    check_det_sequence(data.dets[block_index1], config.get_id());
    check_det_sequence(data.dets[block_index2], config.get_id());
#endif
    LOG("* Rejected \n{}", debug_config_print_end);
  }
} // namespace triqs_cthyb
