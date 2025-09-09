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
    dtau = double(tau_c - tau_c_dag);
    if (histo_proposed) *histo_proposed << dtau;

    // Insert the operators op1 and op2 at time tau1, tau2
    // 1- In the very exceptional case where the insert has failed because an operator is already sitting here
    // (cf std::map doc for insert return), we reject the move.
    // 2- If ok, we store the iterator to the inserted operators for later removal in reject if necessary
    try {
      data.imp_trace.try_insert(tau_c_dag, op_c_dag);
      data.imp_trace.try_insert(tau_c, op_c);
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
    int idx_c_dag, idx_c;
    for (idx_c_dag = 0; idx_c_dag < det_size; ++idx_c_dag) {
      if (det.get_x(idx_c_dag).first < tau_c_dag) break;
    }
    for (idx_c = 0; idx_c < det_size; ++idx_c) {
      if (det.get_y(idx_c).first < tau_c) break;
    }

    // Insert in the det. Returns the ratio of dets (Cf det_manip doc).
    auto det_ratio = det.try_insert(idx_c_dag, idx_c, {tau_c_dag, op_c_dag.inner_index}, {tau_c, op_c.inner_index});

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
    config.insert(tau_c_dag, op_c_dag);
    config.insert(tau_c, op_c);
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
    auto const rs_c_dag = rng(block_size);
    auto const rs_c     = rng(block_size);

    // set operators to be inserted
    op_c_dag = op_desc{
       .block_index = block_index, .inner_index = rs_c_dag, .dagger = true, .linear_index = data.linindex[std::make_pair(block_index, rs_c_dag)]};
    op_c =
       op_desc{.block_index = block_index, .inner_index = rs_c, .dagger = false, .linear_index = data.linindex[std::make_pair(block_index, rs_c)]};

    // choose the tau points uniformly
    tau_c_dag = data.tau_seg.get_random_pt(rng);
    tau_c     = data.tau_seg.get_random_pt(rng);

    // return the proposal probability ratio P_removal / P_insertion
    auto const det_size = static_cast<double>(data.dets[block_index].size());
    return std::pow(block_size * config.beta() / (det_size + 1), 2);
  }

  double move_insert_c_cdag::pauli_proposal() {
    // do a Pauli move with probability pauli_prob
    bool const pauli_move = rng() <= pauli_prob;

    // block determinant and its size
    auto const &det        = data.dets[block_index];
    auto const det_size    = det.size();
    auto const det_size_p1 = static_cast<double>(det_size + 1);

    // choose the inner indices
    auto const rs_c_dag = rng(block_size);
    auto const rs_c     = (pauli_move ? rs_c_dag : rng(block_size));

    // set operators to be inserted
    op_c_dag = op_desc{
       .block_index = block_index, .inner_index = rs_c_dag, .dagger = true, .linear_index = data.linindex[std::make_pair(block_index, rs_c_dag)]};
    op_c =
       op_desc{.block_index = block_index, .inner_index = rs_c, .dagger = false, .linear_index = data.linindex[std::make_pair(block_index, rs_c)]};

    // choose the tau point of c_dag
    tau_c_dag = data.tau_seg.get_random_pt(rng);

    // find operators with the same flavor as c_dag
    c_dag_left_tau.clear(), c_dag_right_tau.clear(), c_left_tau.clear(), c_right_tau.clear();
    for (int i = 0; i < det.size(); ++i) {
      auto const &[tau1, rs1] = det.get_x(i);
      if (rs1 == rs_c_dag) tau1 > tau_c_dag ? c_dag_left_tau.push_back(tau1) : c_dag_right_tau.push_back(tau1);

      auto const &[tau2, rs2] = det.get_y(i);
      if (rs2 == rs_c_dag) tau2 > tau_c_dag ? c_left_tau.push_back(tau2) : c_right_tau.push_back(tau2);
    }

    // move all creation (annihilation) ops to the right of tau1 into c_dag_right_tau (c_right_tau)
    std::ranges::copy(c_dag_left_tau, std::back_inserter(c_dag_right_tau));
    std::ranges::copy(c_left_tau, std::back_inserter(c_right_tau));
    auto const num_c_dag = c_dag_right_tau.size();
    auto const num_c     = c_right_tau.size();

    // insertion/removal proposal probabilities for c_dag
    // - uniform for rs1 and tau1
    // - uniform among all creation operators in the block
    double const p_c_dag_ins = 1.0 / (block_size * config.beta());
    double const p_c_dag_rem = 1.0 / det_size_p1;

    // initialize insertion/removal proposal probabilities for c
    double p_c_ins = 1.0;
    double p_c_rem = 1.0;

    // handle different cases how c can be inserted/removed
    if (rs_c_dag != rs_c) {
      // choose tau point of c uniformly
      tau_c = data.tau_seg.get_random_pt(rng);

      // insertion proposal probability for c:
      // - non-Pauli insertion: uniform for rs2 and tau2
      p_c_ins = (1.0 - pauli_prob) / (block_size * config.beta());

      // removal proposal probability for c:
      // - num_c == 0: uniform among all annihilation operators in the block
      // - num_c > 0: uniform among all annihilation operators in the block only for non-Pauli removals
      p_c_rem = theta(num_c == 0) / det_size_p1 + theta(num_c > 0) * (1.0 - pauli_prob) / det_size_p1;
    } else if (num_c == 0 || num_c_dag == 0) {
      // choose tau point of c uniformly
      tau_c = data.tau_seg.get_random_pt(rng);

      // insertion proposal probability for c:
      // - Pauli insertion: uniform for tau2
      // - non-Pauli insertion: uniform for rs2 and tau2
      p_c_ins = pauli_prob / config.beta() + (1.0 - pauli_prob) / (block_size * config.beta());

      // removal proposal probability for c
      if (num_c == 0) {
        // - Pauli removal: there is only one annihilation operator with the same flavor as c_dag
        // - non-Pauli removal: uniform among all annihilation operators in the block
        p_c_rem = pauli_prob + (1.0 - pauli_prob) / det_size_p1;
      } else {
        // - Pauli removal: c can only be removed if it is a nearest neighbor of c_dag and since num_c > 0, it is chosen
        // uniformly among the left and right nearest neighbor
        // - non-Pauli removal: uniform among all annihilation operators in the block
        bool const pauli_rem = (tau_c_dag - tau_c < tau_c_dag - c_right_tau.front() || tau_c - tau_c_dag < c_right_tau.back() - tau_c_dag);
        p_c_rem              = theta(pauli_rem) * pauli_prob / 2 + (1.0 - pauli_prob) / det_size_p1;
      }
    } else {
      double tau_interval = config.beta();
      bool pauli_ins      = true;
      bool pauli_rem      = true;

      // choose tau point of c uniformly for non-Pauli insertion
      if (!pauli_move) {
        tau_c     = data.tau_seg.get_random_pt(rng);
        pauli_rem = (tau_c_dag - tau_c < tau_c_dag - c_right_tau.front() || tau_c - tau_c_dag < c_right_tau.back() - tau_c_dag);
      }

      // choose tau point for c for Pauli insertion or check if non-Pauli insertion would be a valid Pauli insertion
      if (tau_c_dag - c_dag_right_tau.front() < tau_c_dag - c_right_tau.front()) {
        // the closest rs1 operator to the right of c_dag is a creation operator at tau_right
        auto const tau_right = c_dag_right_tau.front();
        tau_interval         = static_cast<double>(tau_c_dag - tau_right);

        if (pauli_move) {
          // choose tau point of c uniformly between tau_right and tau1
          tau_c = tau_right + data.tau_seg.get_random_pt(rng, tau_c_dag - tau_right);
        } else {
          // is the non-Pauli insertion a valid Pauli insertion?
          pauli_ins = (tau_c_dag - tau_c < tau_c_dag - tau_right);
        }
      } else {
        // the closest rs1 operator to the left of c_dag is at tau_left (creation or annihilation)
        auto const tau_left = (c_dag_right_tau.back() - tau_c_dag < c_right_tau.back() - tau_c_dag ? c_dag_right_tau.back() : c_right_tau.back());
        tau_interval        = static_cast<double>(tau_left - tau_c_dag);

        if (pauli_move) {
          // choose tau point of c uniformly between tau1 and tau_left
          tau_c = tau_c_dag + data.tau_seg.get_random_pt(rng, tau_left - tau_c_dag);
        } else {
          // is the non-Pauli insertion a valid Pauli insertion?
          pauli_ins = (tau_c - tau_c_dag < tau_left - tau_c_dag);
        }
      }

      // insertion proposal probability for c:
      // - Pauli insertion: uniform for tau2 on the tau_interval (only if it a valid Pauli insertion)
      // - non-Pauli insertion: uniform for rs2 and tau2
      p_c_ins = theta(pauli_ins) * pauli_prob / tau_interval + (1.0 - pauli_prob) / (block_size * config.beta());

      // removal proposal probability for c:
      // - Pauli removal: c is a nearest neighbor of c_dag and num_c > 0, so it is chosen uniformly among the left and
      // right nearest neighbor
      // - non-Pauli removal: uniform among all annihilation operators in the block
      p_c_rem = theta(pauli_rem) * pauli_prob / 2 + (1. - pauli_prob) / det_size_p1;
    }

    return (p_c_dag_rem * p_c_rem) / (p_c_dag_ins * p_c_ins);
  }

} // namespace triqs_cthyb
