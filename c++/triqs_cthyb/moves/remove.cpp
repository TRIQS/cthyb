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

#include "./remove.hpp"
#include "./util.hpp"

namespace triqs_cthyb {

  move_remove_c_cdag::move_remove_c_cdag(int block_index, int block_size, std::string const &block_name, qmc_data &data,
                                         mc_tools::random_generator &rng, histo_map_t *histos)
     : data(data),
       config(data.config),
       rng(rng),
       block_index(block_index),
       block_size(block_size),
       histo_proposed(add_histo("remove_length_proposed_" + block_name, histos, config.beta())),
       histo_accepted(add_histo("remove_length_accepted_" + block_name, histos, config.beta())) {}

  // -----------------------------------------------------------

  mc_weight_t move_remove_c_cdag::attempt() {

    LOG("{}\n* Attempt for remove_c_cdag (block {})", debug_config_print_start, block_index);

    auto &det = data.dets[block_index];

    // Pick up a couple of C, Cdagger to remove at random
    // Remove the operators from the traces
    int det_size = det.size();
    if (det_size == 0) return 0; // nothing to remove
    int num_c_dag = rng(det_size), num_c = rng(det_size);
    LOG("* Proposing to remove for block {}:\n       {}-th Cdag\n       {}-th C", block_index, num_c_dag, num_c);

    // now mark 2 nodes for deletion
    // FIXME rename tau1 -> tau_c etc...
    tau1 = det.get_y(num_c).first;
    tau2 = det.get_x(num_c_dag).first;

    data.imp_trace.try_delete(tau1);
    data.imp_trace.try_delete(tau2);

    // record the length of the proposed removal
    dtau = double(tau2 - tau1);
    if (histo_proposed) *histo_proposed << dtau;

    auto det_ratio = det.try_remove(num_c_dag, num_c);

    // proposition probability
    auto t_ratio = std::pow(block_size * config.beta() / double(det_size), 2); // Size of the det before the try_delete!

    // For quick abandon
    double random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double p_yee = std::abs(det_ratio / t_ratio / data.atomic_weight);

    // recompute the atomic_weight
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) {
      LOG("atomic_weight == 0");
      return 0;
    }
    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(remove) atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                          << new_atomic_weight / data.atomic_weight << " in config " << config.get_id();

    mc_weight_t p = atomic_weight_ratio * det_ratio;

    LOG("Trace ratio: {}  Det ratio: {}  Prefactor: {} ", atomic_weight_ratio, det_ratio, t_ratio);
    LOG("Weight: {}\n  p_yee * newtrace: {} ", p / t_ratio, p_yee * new_atomic_weight);
    ALWAYS_EXPECTS(isfinite(p), "(remove) p = {} not finite", p); 
    ALWAYS_EXPECTS(isfinite(p), "(remove) p / t_ratio  not finite. p = {}, t_ratio = {} ", p, t_ratio); 

    return p / t_ratio;
  }

  // -----------------------------------------------------------

  mc_weight_t move_remove_c_cdag::accept() {

    // remove from the tree
    data.imp_trace.confirm_delete();

    // remove from the configuration
    config.erase(tau1);
    config.erase(tau2);

    // remove from the determinants
    data.dets[block_index].complete_operation();
    data.update_sign();
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;
    if (histo_accepted) *histo_accepted << dtau;

#ifdef EXT_DEBUG
    check_det_sequence(data.dets[block_index], config.get_id());
#endif

    LOG("* Accepted \n{}", debug_config_print_end);
    return data.current_sign / data.old_sign;
  }

  // -----------------------------------------------------------

  void move_remove_c_cdag::reject() {

    data.imp_trace.cancel_delete();
    data.dets[block_index].reject_last_try();

#ifdef EXT_DEBUG
    check_det_sequence(data.dets[block_index], config.get_id());
#endif
    LOG("* Rejected \n{}", debug_config_print_end);
  }
} // namespace triqs_cthyb
