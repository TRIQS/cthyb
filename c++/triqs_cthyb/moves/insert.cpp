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

namespace triqs_cthyb {

  histogram *move_insert_c_cdag::add_histo(std::string const &name, histo_map_t *histos) {
    if (!histos) return nullptr;
    auto new_histo = histos->insert({name, {.0, config.beta(), 100}});
    return &(new_histo.first->second);
  }

  move_insert_c_cdag::move_insert_c_cdag(int block_index, int block_size, std::string const &block_name, qmc_data &data,
                                         mc_tools::random_generator &rng, histo_map_t *histos, double pauli_prob)
     : data(data),
       config(data.config),
       rng(rng),
       block_index(block_index),
       block_size(block_size),
       histo_proposed(add_histo("insert_length_proposed_" + block_name, histos)),
       histo_accepted(add_histo("insert_length_accepted_" + block_name, histos)),
       pauli_prob(pauli_prob) {}

  mc_weight_t move_insert_c_cdag::attempt() {

#ifdef EXT_DEBUG
    std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cerr << "In config " << config.get_id() << std::endl;
    std::cerr << "* Attempt for move_insert_c_cdag (block " << block_index << ")" << std::endl;
#endif

    // propose tau points and inner block indices --> set the proposal distribution ratio
    double t_ratio = (pauli_prob <= 0.0 ? uniform_proposal() : pauli_proposal());

#ifdef EXT_DEBUG
    std::cerr << "* Proposing to insert:" << std::endl;
    std::cerr << op1 << " at " << tau1 << std::endl;
    std::cerr << op2 << " at " << tau2 << std::endl;
#endif

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
      std::cerr << "Insert error : recovering ... " << std::endl;
      data.imp_trace.cancel_insert();
      return 0;
    }

    // Computation of det ratio
    auto &det    = data.dets[block_index];
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
      TRIQS_RUNTIME_ERROR << "(insert) trace_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                          << new_atomic_weight / data.atomic_weight << " in config " << config.get_id();

    mc_weight_t p = atomic_weight_ratio * det_ratio;

#ifdef EXT_DEBUG
    std::cerr << "Atomic ratio: " << atomic_weight_ratio << '\t';
    std::cerr << "Det ratio: " << det_ratio << '\t';
    std::cerr << "Prefactor: " << t_ratio << '\t';
    std::cerr << "Weight: " << p * t_ratio << std::endl;
    std::cerr << "p_yee * newtrace: " << p_yee * new_atomic_weight << std::endl;
#endif

    if (!isfinite(p * t_ratio)) {
      std::cerr << "Insert move info:\n";
      std::cerr << "Atomic ratio: " << atomic_weight_ratio << '\t';
      std::cerr << "Det ratio: " << det_ratio << '\t';
      std::cerr << "Prefactor: " << t_ratio << '\t';
      std::cerr << "Weight: " << p * t_ratio << std::endl;
      std::cerr << "p_yee * newtrace: " << p_yee * new_atomic_weight << std::endl;

      TRIQS_RUNTIME_ERROR << "(insert) p * t_ratio not finite p : " << p << " t_ratio : " << t_ratio << " in config " << config.get_id();
    }
    return p * t_ratio;
  }

  mc_weight_t move_insert_c_cdag::accept() {

    // insert in the tree
    data.imp_trace.confirm_insert();

    // insert in the configuration
    config.insert(tau1, op1);
    config.insert(tau2, op2);
    config.finalize();

    // insert in the determinant
    data.dets[block_index].complete_operation();
    data.update_sign();
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;
    if (histo_accepted) *histo_accepted << dtau;

#ifdef EXT_DEBUG
    std::cerr << "* Move move_insert_c_cdag accepted" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    check_det_sequence(data.dets[block_index], config.get_id());
#endif

    return data.current_sign / data.old_sign;
  }

  void move_insert_c_cdag::reject() {
    config.finalize();
    data.imp_trace.cancel_insert();
    data.dets[block_index].reject_last_try();

#ifdef EXT_DEBUG
    std::cerr << "* Move move_insert_c_cdag rejected" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    check_det_sequence(data.dets[block_index], config.get_id());
#endif
  }

  double move_insert_c_cdag::uniform_proposal() {
    // choose the inner block indices uniformly
    auto const rs1 = rng(block_size);
    auto const rs2 = rng(block_size);
    op1 = op_desc{.block_index = block_index, .inner_index = rs1, .dagger = true, .linear_index = data.linindex[std::make_pair(block_index, rs1)]};
    op2 = op_desc{.block_index = block_index, .inner_index = rs2, .dagger = false, .linear_index = data.linindex[std::make_pair(block_index, rs2)]};

    // choose the tau points uniformly
    tau1 = data.tau_seg.get_random_pt(rng);
    tau2 = data.tau_seg.get_random_pt(rng);

    // return the proposal probability ratio P_removal / P_insertion
    auto const det_size = static_cast<double>(data.dets[block_index].size());
    return std::pow(block_size * config.beta() / (det_size + 1), 2);
  }

  double move_insert_c_cdag::pauli_proposal() {
    // do a Pauli move with probability pauli_prob
    bool const pauli_move = rng() <= pauli_prob;

    // choose the inner indices and initialize operators to be inserted
    auto const rs1 = rng(block_size);
    auto const rs2 = (pauli_move ? rs1 : rng(block_size));
    op1 = op_desc{.block_index = block_index, .inner_index = rs1, .dagger = true, .linear_index = data.linindex[std::make_pair(block_index, rs1)]};
    op2 = op_desc{.block_index = block_index, .inner_index = rs2, .dagger = false, .linear_index = data.linindex[std::make_pair(block_index, rs2)]};

    // block determinant and its size
    auto const &det        = data.dets[block_index];
    auto const det_size    = det.size();
    auto const det_size_p1 = static_cast<double>(det_size + 1);

    // choose the tau point of c_dag and c (only for non-Pauli moves)
    tau1 = data.tau_seg.get_random_pt(rng);
    if (!pauli_move) tau2 = data.tau_seg.get_random_pt(rng);

    // initialize the proposal probability ratio P_removal / P_insertion
    double t_ratio = block_size * config.beta() / det_size_p1;

    // find operators with the same flavor as c_dag
    c_dag_left_tau.clear(), c_dag_right_tau.clear(), c_left_tau.clear(), c_right_tau.clear();
    for (int i = 0; i < det.size(); ++i) {
      auto const &[tau_dag, rs_dag] = det.get_x(i);
      if (rs_dag == rs1) tau_dag > tau1 ? c_dag_left_tau.push_back(tau_dag) : c_dag_right_tau.push_back(tau_dag);

      auto const &[tau, rs] = det.get_y(i);
      if (rs == rs1) tau > tau1 ? c_left_tau.push_back(tau) : c_right_tau.push_back(tau);
    }

    // move all creation (annihilation) ops to the right of tau1 into c_dag_right_tau (c_right_tau)
    std::ranges::copy(c_dag_left_tau, std::back_inserter(c_dag_right_tau));
    std::ranges::copy(c_left_tau, std::back_inserter(c_right_tau));

    // choose tau point for c and update P_insertion of t_ratio
    if (!c_dag_right_tau.empty() && !c_right_tau.empty() && rs1 == rs2) {
      // rs1 == rs2 && at least one creation and annihilation operator with rs1 is already present
      if (tau1 - c_dag_right_tau.front() < tau1 - c_right_tau.front()) {
        // the closest rs1-op to the right of tau1 is a creation operator at tau_right
        auto const tau_right = c_dag_right_tau.front();

        // for Pauli moves, choose tau2 between tau_right and tau1
        if (pauli_move) tau2 = tau_right + data.tau_seg.get_random_pt(rng, tau1 - tau_right);

        // update t_ratio
        if ((tau1 - tau2) < (tau1 - tau_right)) {
          // c is inserted between c_dag and the closest rs1-creation op to the right of tau1 (can be Pauli move or non-Pauli move)
          t_ratio /= pauli_prob / static_cast<double>(tau1 - tau_right) + (1. - pauli_prob) / (block_size * config.beta());
        } else {
          // c is inserted somewhere else (can only be a non-Pauli move)
          t_ratio *= block_size * config.beta() / (1. - pauli_prob);
        }
      } else {
        // the closest rs1-op to the left of tau1 is at tau_left
        auto const tau_left = (c_dag_right_tau.back() - tau1 < c_right_tau.back() - tau1 ? c_dag_right_tau.back() : c_right_tau.back());

        // for Pauli moves, choose tau2 between tau1 and tau_left
        if (pauli_move) tau2 = tau1 + data.tau_seg.get_random_pt(rng, tau_left - tau1);

        // update t_ratio
        if ((tau2 - tau1) < (tau_left - tau1)) {
          // c is inserted between c_dag and the closest rs1-op to the left of tau1 (can be Pauli move or non-Pauli move)
          t_ratio /= pauli_prob / static_cast<double>(tau_left - tau1) + (1. - pauli_prob) / (block_size * config.beta());
        } else {
          // c is inserted somewhere else (can only be a non-Pauli move)
          t_ratio *= block_size * config.beta() / (1. - pauli_prob);
        }
      }
    } else {
      // rs1 != rs2 || no creation or no annihilation operator with rs1 is present
      // for Pauli moves, choose tau2 uniformly on [0, beta)
      if (pauli_move) tau2 = data.tau_seg.get_random_pt(rng);

      // update t_ratio
      if (rs1 == rs2) {
        // can be Pauli move or non-Pauli move
        t_ratio /= (pauli_prob + (1. - pauli_prob) / double(block_size)) / config.beta();
      } else {
        // can only be a non-Pauli move
        t_ratio *= block_size * config.beta() / (1. - pauli_prob);
      }
    }

    // update P_removal of t_ratio
    auto const num_pauli = std::min(c_right_tau.size() + (rs1 == rs2 ? 1 : 0), 2ul);
    if (num_pauli == 0 || num_pauli == det_size + 1) {
      // num_pauli == 0 only if rs1 != rs2 (only non-Pauli moves), num_pauli == det_size + 1 only if rs1 == rs2 (Pauli and non-Pauli moves)
      t_ratio /= det_size_p1;
    } else {
      // find the closest rs1-annihilation operators to tau1 (can be tau2 or some existing operator)
      auto tau_right = (c_right_tau.empty() ? tau2 : c_right_tau.front());
      auto tau_left  = (c_right_tau.empty() ? tau2 : c_right_tau.back());
      if (!c_right_tau.empty()) {
        if ((tau1 - tau2) < (tau1 - tau_right) && (rs1 == rs2)) tau_right = tau2;
        if ((tau2 - tau1) < (tau_left - tau1) && (rs1 == rs2)) tau_left = tau2;
      }
      if (tau2 == tau_right || tau2 == tau_left) {
        // tau2 is the right/left nearest neighbor of tau1 (can be Pauli move or non-Pauli move)
        t_ratio *= pauli_prob / double(num_pauli) + (1. - pauli_prob) / det_size_p1;
      } else {
        // tau2 is not the nearest neighbor of tau1 (can only be a non-Pauli move)
        t_ratio *= (1. - pauli_prob) / det_size_p1;
      }
    }

    return t_ratio;
  }

} // namespace triqs_cthyb
