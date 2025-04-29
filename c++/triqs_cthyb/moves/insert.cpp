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
       pauli_prob(pauli_prob),
       Nmax(100) {
         if (pauli_prob > 0.0) vec_ind.reserve(Nmax);
       }

  mc_weight_t move_insert_c_cdag::attempt() {

#ifdef EXT_DEBUG
    std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cerr << "In config " << config.get_id() << std::endl;
    std::cerr << "* Attempt for move_insert_c_cdag (block " << block_index << ")" << std::endl;
#endif

    bool pauli_move = false;
    if (pauli_prob > 0.0) {
      double ran = rng();
      if (ran <= pauli_prob) pauli_move = true;
    }

    // Pick up the value of alpha and choose the operators
    auto rs1 = rng(block_size);
    auto rs2 = (pauli_move ? rs1 : rng(block_size));
    op1 = op_desc{block_index, rs1, true, data.linindex[std::make_pair(block_index, rs1)]};
    op2 = op_desc{block_index, rs2, false, data.linindex[std::make_pair(block_index, rs2)]};

    auto &det    = data.dets[block_index];
    int det_size = det.size();

    // Choice of times for insertion. Find the time as double and them put them on the grid.
    tau1 = data.tau_seg.get_random_pt(rng);
    if (!pauli_move) tau2 = data.tau_seg.get_random_pt(rng);
    mc_weight_t t_ratio = 1.;

    if (pauli_prob > 0.0) {

      t_ratio = block_size * config.beta() / double(det_size + 1);

      int ic_dagL = -1, ic_dagR = -1, ic_nodagL = -1, ic_nodagR = -1;

      if (det_size > Nmax) {
        while (det_size > Nmax) Nmax *= 2;
        vec_ind.reserve(Nmax);
      }

      vec_ind.clear();
      int j = 0, op_pos = -1;

      if (rs1 == rs2) {
        // Look at position of tau1 among creation operators of the same flavor
        for (int i = 0; i < det_size; ++i) {
          if (det.get_x(i).first < tau1 && op_pos == -1) op_pos = j;
          if (det.get_x(i).second != rs1) continue;
          vec_ind.push_back(i);
          ++j;
        }
      }

      int size = vec_ind.size();

      if (op_pos == -1) op_pos = size;

      // Creation operators at the left and right of tau1
      if (size != 0) {
        ic_dagR = (op_pos == size ? vec_ind[0] : vec_ind[op_pos]);
        ic_dagL = (op_pos == 0 ? vec_ind[size-1] : vec_ind[op_pos-1]);
      }

      vec_ind.clear();
      j = 0, op_pos = -1;

      // Look at position of tau1 among annihilation operators of the same flavor
      for (int i = 0; i < det_size; ++i) {
        if (det.get_y(i).first < tau1 && op_pos == -1) op_pos = j;
        if (det.get_y(i).second != rs1) continue;
        vec_ind.push_back(i);
        ++j;
      }

      size = vec_ind.size();

      if (op_pos == -1) op_pos = size;

      time_pt tRnodag, tLnodag;

      // Annihilation operators at the left and right of tau1
      if (size != 0) {
        ic_nodagR = (op_pos == size ? vec_ind[0] : vec_ind[op_pos]);
        ic_nodagL = (op_pos == 0 ? vec_ind[size-1] : vec_ind[op_pos-1]);
        tRnodag = det.get_y(ic_nodagR).first;
        tLnodag = det.get_y(ic_nodagL).first;
      }

      if (ic_nodagR != -1 && ic_dagR != -1) {

        auto tRdag = det.get_x(ic_dagR).first;
        auto tLdag = det.get_x(ic_dagL).first;

        auto tR = ((tau1 - tRdag) > (tau1 - tRnodag) ? tRnodag : tRdag);
        auto tL = ((tLdag - tau1) > (tLnodag - tau1) ? tLnodag : tLdag);

        if (tR == tRdag) {
          if (pauli_move) tau2 = tR + data.tau_seg.get_random_pt(rng, tau1 - tR);
          if ((tau1 - tau2) < (tau1 - tR))
            t_ratio /= pauli_prob / double(tau1 - tR) + (1. - pauli_prob) / (block_size * config.beta());
          else
            t_ratio *= block_size * config.beta() / (1. - pauli_prob);
        }
        else {
          if (pauli_move) tau2 = tau1 + data.tau_seg.get_random_pt(rng, tL - tau1);
          if ((tau2 - tau1) < (tL - tau1))
            t_ratio /= pauli_prob / double(tL - tau1) + (1. - pauli_prob) / (block_size * config.beta());
          else
            t_ratio *= block_size * config.beta() / (1. - pauli_prob);
        }
      }
      else { // if no operators of the same flavor or different flavors for insertion
        if (pauli_move) tau2 = data.tau_seg.get_random_pt(rng);
        if (rs1 == rs2)
          t_ratio /= (pauli_prob + (1. - pauli_prob) / double(block_size)) / config.beta();
        else
          t_ratio *= block_size * config.beta() / (1. - pauli_prob);
      }

      int num_pauli = size;
      if (rs1 == rs2) ++num_pauli;
      num_pauli = std::min(num_pauli, 2);
      if (num_pauli == 0 || num_pauli == det_size + 1)
        t_ratio /= double(det_size + 1);
      else {
        if (size == 0) {   // In this case rs1 = rs2
          tRnodag = tau2;
          tLnodag = tau2;
        }
        else {
          if ((tau1 - tau2) < (tau1 - tRnodag) && (rs1 == rs2)) tRnodag = tau2;
          if ((tau2 - tau1) < (tLnodag - tau1) && (rs1 == rs2)) tLnodag = tau2;
        }
        if (tau2 == tRnodag || tau2 == tLnodag)
          t_ratio *= pauli_prob / double(num_pauli) + (1. - pauli_prob) / double(det_size + 1);
        else
          t_ratio *= (1. - pauli_prob) / double(det_size + 1);
      }
    }

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
    if (pauli_prob == 0.0) t_ratio = std::pow(block_size * config.beta() / double(det.size() + 1), 2);

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
}
