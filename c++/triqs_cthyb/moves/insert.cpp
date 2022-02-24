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

#include "./insert.hpp"
#include "./util.hpp"

namespace triqs_cthyb {

  move_insert_c_cdag::move_insert_c_cdag(int block_index, int block_size, std::string const &block_name, qmc_data &data,
                                         mc_tools::random_generator &rng, histo_map_t *histos)
     : data(data),
       config(data.config),
       rng(rng),
       block_index(block_index),
       block_size(block_size),
       histo_proposed(add_histo("insert_length_proposed_" + block_name, histos, config.beta())),
       histo_accepted(add_histo("insert_length_accepted_" + block_name, histos, config.beta())) {}

  // -----------------------------------------------------------
  
  mc_weight_t move_insert_c_cdag::attempt() {

    SPDLOG_TRACE("{}\n* Attempt for move_insert_c_cdag (block {})", debug_config_print_start, block_index);

    // Pick up the value of alpha and choose the operators
    auto rs1 = rng(block_size), rs2 = rng(block_size);
    op1 = op_desc{block_index, rs1, true, data.linindex[std::make_pair(block_index, rs1)]};
    op2 = op_desc{block_index, rs2, false, data.linindex[std::make_pair(block_index, rs2)]};

    // Choice of times for insertion. Find the time as double and them put them on the grid.
    tau1 = data.tau_seg.get_random_pt(rng);
    tau2 = data.tau_seg.get_random_pt(rng);

    LOG("* Proposing to insert:\n   {} at {}\n   {} at {}", op1, tau1, op2, tau2);

    // record the length of the proposed insertion
    dtau = double(tau2 - tau1);
    if (histo_proposed) *histo_proposed << dtau;

    // Insert the operators op1 and op2 at time tau1, tau2
    // 1- In the very exceptional case where the insert has failed because an operator is already sitting here
    // (cf std::map doc for insert return), we reject the move.
    // 2- If ok, we store the iterator to the inserted operators for later removal in reject if necessary
    try {
      data.imp_trace.try_insert(tau1, op1);
      data.imp_trace.try_insert(tau2, op2);
    } catch (rbt_insert_error const &) {
      LOG( "Insert error : recovering ... ");
      data.imp_trace.cancel_insert();
      return 0;
    }

    // Computation of det ratio
    auto &det    = data.dets[block_index];
    int det_size = det.size();

    // Find the position for insertion in the determinant
    // NB : the determinant stores the C in decreasing time order.
    // FIXME : not great ... should be binary search
    // but hte detmanip doesnot have the list in order

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
      LOG("atomic_weight == 0");
      return 0;
    }
    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(insert) trace_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                          << new_atomic_weight / data.atomic_weight << " in config " << config;

    mc_weight_t p = atomic_weight_ratio * det_ratio;

    LOG("Atomic ratio: {}  Det ratio: {}  Prefactor: {} ", atomic_weight_ratio, det_ratio, t_ratio);
    LOG("Weight: {}\n  p_yee * newtrace: {} ", p * t_ratio, p_yee * new_atomic_weight);
    ALWAYS_EXPECTS(isfinite(p * t_ratio), "(insert) p * t_ratio not finite p = {}, t_ratio = {} \n config:\n ", p, t_ratio, config);
    
    return p * t_ratio;
  }

  // -----------------------------------------------------------

  mc_weight_t move_insert_c_cdag::accept() {

    // insert in the tree
    data.imp_trace.confirm_insert();

    // insert in the configuration
    config.insert(tau1, op1);
    config.insert(tau2, op2);

    // insert in the determinant
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

  void move_insert_c_cdag::reject() {

    data.imp_trace.cancel_insert();
    data.dets[block_index].reject_last_try();

#ifdef EXT_DEBUG
    check_det_sequence(data.dets[block_index], config.get_id());
#endif
    LOG("* Rejected \n{}", debug_config_print_end);
  }
} // namespace triqs_cthyb
