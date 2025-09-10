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
#include "../config.hpp"
#include "../qmc_data.hpp"
#include "../types.hpp"

#include <triqs/mc_tools.hpp>
#include <triqs/stat/histograms.hpp>
#include <triqs/utility/exceptions.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <string>
#include <tuple>

namespace triqs_cthyb {

  histogram *move_remove_c_cdag::add_histo(std::string const &name, histo_map_t *histos) {
    if (!histos) return nullptr;
    auto new_histo = histos->insert({name, {.0, config.beta(), 100}});
    return &(new_histo.first->second);
  }

  move_remove_c_cdag::move_remove_c_cdag(int block_index, int block_size, std::string const &block_name, qmc_data &data,
                                         mc_tools::random_generator &rng, histo_map_t *histos, double pauli_prob)
     : data(data),
       config(data.config),
       rng(rng),
       block_index(block_index),
       block_size(block_size),
       histo_proposed(add_histo("remove_length_proposed_" + block_name, histos)),
       histo_accepted(add_histo("remove_length_accepted_" + block_name, histos)),
       pauli_prob(pauli_prob) {}

  mc_weight_t move_remove_c_cdag::attempt() {

#ifdef EXT_DEBUG
    std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cerr << "In config " << config.get_id() << std::endl;
    std::cerr << "* Attempt for move_remove_c_cdag (block " << block_index << ")" << std::endl;
#endif

    // block determinant and its size
    auto &det           = data.dets[block_index];
    auto const det_size = static_cast<int>(det.size());

    // early return if nothing to remove
    if (det_size == 0) return 0;

    // propose operators to remove and set the proposal distribution ratio
    double const t_ratio = (pauli_prob <= 0.0 ? uniform_proposal() : pauli_proposal());

#ifdef EXT_DEBUG
    std::cerr << "* Proposing to remove: ";
    std::cerr << num_c_dag << "-th Cdag(" << block_index << ",...), ";
    std::cerr << num_c << "-th C(" << block_index << ",...)" << std::endl;
#endif

    // mark the operators for deletion in the impurity trace
    tau_c     = data.imp_trace.try_delete(idx_c, block_index, false);
    tau_c_dag = data.imp_trace.try_delete(idx_c_dag, block_index, true);

    // gather performance statistics - record the length of the proposed removal
    dtau = double(tau_c_dag - tau_c);
    if (histo_proposed) *histo_proposed << dtau;

    // remove the ops from the determinant and get the determinant ratio
    auto const det_ratio = det.try_remove(idx_c_dag, idx_c);

    // for early rejection
    double const random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double const p_yee = std::abs(det_ratio / t_ratio / data.atomic_weight);

    // computation of the new/old impurity trace
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) {
#ifdef EXT_DEBUG
      std::cerr << "atomic_weight == 0" << std::endl;
#endif
      return 0;
    }

    // impurity trace ratio
    auto const atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(remove) atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                          << new_atomic_weight / data.atomic_weight << " in config " << config.get_id();

    // weight ratio
    mc_weight_t const p = atomic_weight_ratio * det_ratio;

#ifdef EXT_DEBUG
    std::cerr << "Trace ratio: " << atomic_weight_ratio << '\t';
    std::cerr << "Det ratio: " << det_ratio << '\t';
    std::cerr << "Prefactor: " << t_ratio << '\t';
    std::cerr << "Weight: " << p / t_ratio << std::endl;
#endif

    if (!isfinite(p)) {
      std::cerr << "Remove move info\n";
      std::cerr << "Trace ratio: " << atomic_weight_ratio << '\t';
      std::cerr << "Det ratio: " << det_ratio << '\t';
      std::cerr << "Prefactor: " << t_ratio << '\t';
      std::cerr << "Weight: " << p / t_ratio << std::endl;
      TRIQS_RUNTIME_ERROR << "(remove) p not finite :" << p << " in config " << config.get_id();
    }

    if (!isfinite(p / t_ratio)) {
      TRIQS_RUNTIME_ERROR << "(remove) p / t_ratio not finite p : " << p << " t_ratio :  " << t_ratio << " in config " << config.get_id();
    }

    return p / t_ratio;
  }

  mc_weight_t move_remove_c_cdag::accept() {

    // confirm the removal from the impurity trace
    data.imp_trace.confirm_delete();

    // remove from the configuration
    config.erase(tau_c);
    config.erase(tau_c_dag);
    config.finalize();

    // complete the removal from the determinant
    data.dets[block_index].complete_operation();
    data.update_sign();
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;

    // gather performance statistics
    if (histo_accepted) *histo_accepted << dtau;

#ifdef EXT_DEBUG
    std::cerr << "* Move move_remove_c_cdag accepted" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    check_det_sequence(data.dets[block_index], config.get_id());
#endif

    // return sign correction
    return static_cast<double>(data.current_sign) / data.old_sign;
  }

  void move_remove_c_cdag::reject() {
    // reject insertions into the impurity trace and determinant
    config.finalize();
    data.imp_trace.cancel_delete();
    data.dets[block_index].reject_last_try();

#ifdef EXT_DEBUG
    std::cerr << "* Move move_remove_c_cdag rejected" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    check_det_sequence(data.dets[block_index], config.get_id());
#endif
  }

  double move_remove_c_cdag::uniform_proposal() {
    // block determinant and its size
    auto const &det     = data.dets[block_index];
    auto const det_size = static_cast<int>(det.size());

    // choose random creation and annihilation operators to remove
    idx_c_dag = rng(det_size);
    idx_c     = rng(det_size);

    // return the proposal probability ratio P_removal / P_insertion
    return std::pow(block_size * config.beta() / static_cast<double>(det_size), 2);
  }

  double move_remove_c_cdag::pauli_proposal() {
    // block determinant and its size
    auto const &det     = data.dets[block_index];
    auto const det_size = static_cast<int>(det.size());

    // choose a random c_dag to remove
    idx_c_dag           = rng(det_size);
    tau_c_dag           = det.get_x(idx_c_dag).first;
    auto const rs_c_dag = det.get_x(idx_c_dag).second;

    // find operators with the same flavor as c_dag (c_dag is c_dag_right[0])
    c_dag_left.clear(), c_dag_right.clear(), c_left.clear(), c_right.clear();
    for (int i = 0; i < det_size; ++i) {
      auto const &[tau1, rs1] = det.get_x(i);
      if (rs1 == rs_c_dag) tau1 > tau_c_dag ? c_dag_left.push_back(i) : c_dag_right.push_back(i);

      auto const &[tau2, rs2] = det.get_y(i);
      if (rs2 == rs_c_dag) tau2 > tau_c_dag ? c_left.push_back(i) : c_right.push_back(i);
    }

    // move all creation (annihilation) ops to the right of c_dag into c_dag_right (c_right)
    std::ranges::copy(c_dag_left, std::back_inserter(c_dag_right));
    std::ranges::copy(c_left, std::back_inserter(c_right));
    auto const num_c_dag = c_dag_right.size();
    auto const num_c     = c_right.size();
    auto const num_pauli = std::min(2ul, num_c);

    // insertion/removal proposal probabilities for c_dag (see insert move)
    double const p_c_dag_ins = 1.0 / (block_size * config.beta());
    double const p_c_dag_rem = 1.0 / static_cast<double>(det_size);

    // initialize insertion/removal proposal probabilities for c
    double p_c_ins = 1.0;
    double p_c_rem = 1.0;

    // choose a c to remove and set the removal proposal probability
    if (num_pauli > 0) {
      // do a Pauli move with probability pauli_prob
      bool const pauli_move = (num_pauli == det_size || rng() <= pauli_prob);

      // choose c and check if the chosen c be removed with a Pauli removal (always true for Pauli moves)
      bool pauli_rem = true;
      if (pauli_move) {
        // for Pauli removals, choose between the c to the left and right of c_dag
        idx_c = (num_pauli == 1 ? c_right.front() : (rng(2) == 0 ? c_right.front() : c_right.back()));
      } else {
        // for non-Pauli removals, choose c uniformly among all annihilation operators in the block
        idx_c     = rng(det_size);
        pauli_rem = (idx_c == c_right.front() || idx_c == c_right.back());
      }
      p_c_rem = theta(pauli_rem) * pauli_prob / static_cast<double>(num_pauli) + (1. - pauli_prob) / static_cast<double>(det_size);
    } else {
      // if there is no annihilation operator of the same flavor as c_dag, choose c uniformly
      idx_c   = rng(det_size);
      p_c_rem = 1.0 / static_cast<double>(det_size);
    }
    tau_c           = det.get_y(idx_c).first;
    auto const rs_c = det.get_y(idx_c).second;

    // set the insertion proposal probability
    if (rs_c != rs_c_dag) {
      // c_dag and c are of different flavor (only non-Pauli insertions)
      p_c_ins = (1. - pauli_prob) / (block_size * config.beta());
    } else if (num_c == 1 || num_c_dag == 1) {
      // there is no other c_dag or c of the same flavor as c_dag
      p_c_ins = pauli_prob / config.beta() + (1. - pauli_prob) / (block_size * config.beta());
    } else {
      // tau points of right/left nearest neighbors of c_dag of the same flavor (excluding c_dag and c)
      auto const tau_c_right     = (idx_c == c_right[0] ? det.get_y(c_right[1]).first : det.get_y(c_right[0]).first);
      auto const tau_c_left      = (idx_c == c_right[num_c - 1] ? det.get_y(c_right[num_c - 2]).first : det.get_y(c_right[num_c - 1]).first);
      auto const tau_c_dag_right = det.get_x(c_dag_right[1]).first;
      auto const tau_c_dag_left  = det.get_x(c_dag_right[num_c_dag - 1]).first;

      // check if c could have been inserted with a Pauli insertion and get relevant tau interval (see insert move)
      bool pauli_ins      = false;
      double tau_interval = config.beta();
      if (tau_c_dag - tau_c_dag_right < tau_c_dag - tau_c_right) {
        // consider tau interval to the right of c_dag
        tau_interval = static_cast<double>(tau_c_dag - tau_c_dag_right);
        pauli_ins    = (tau_c_dag - tau_c < tau_c_dag - tau_c_dag_right);
      } else {
        // consider tau interval to the left of c_dag
        auto const tau_left = (tau_c_dag_left - tau_c_dag < tau_c_left - tau_c_dag ? tau_c_dag_left : tau_c_left);
        tau_interval        = static_cast<double>(tau_left - tau_c_dag);
        pauli_ins           = (tau_c - tau_c_dag < tau_left - tau_c_dag);
      }

      p_c_ins = theta(pauli_ins) * pauli_prob / tau_interval + (1. - pauli_prob) / (block_size * config.beta());
    }

    return (p_c_dag_rem * p_c_rem) / (p_c_dag_ins * p_c_ins);
  }

} // namespace triqs_cthyb
