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

namespace triqs_cthyb {

  histogram * move_remove_c_cdag::add_histo(std::string const &name, histo_map_t *histos) {
    if (!histos) return nullptr;
    auto new_histo = histos->insert({name, {.0, config.beta(), 100}});
    return &(new_histo.first->second);
  }

  move_remove_c_cdag::move_remove_c_cdag(int block_index, int block_size, std::string const &block_name, qmc_data &data, mc_tools::random_generator &rng,
                                         histo_map_t *histos, double pauli_prob)
     : data(data),
       config(data.config),
       rng(rng),
       block_index(block_index),
       block_size(block_size),
       histo_proposed(add_histo("remove_length_proposed_" + block_name, histos)),
       histo_accepted(add_histo("remove_length_accepted_" + block_name, histos)),
       Nmax(100),
       pauli_prob(pauli_prob) {
         if (pauli_prob > 0.0) vec_ind.reserve(Nmax);
       }

  mc_weight_t move_remove_c_cdag::attempt() {

#ifdef EXT_DEBUG
    std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cerr << "In config " << config.get_id() << std::endl;
    std::cerr << "* Attempt for move_remove_c_cdag (block " << block_index << ")" << std::endl;
#endif

    auto &det = data.dets[block_index];

    // Pick up a couple of C, Cdagger to remove at random
    // Remove the operators from the traces
    int det_size = det.size();
    if (det_size == 0) return 0; // nothing to remove
    int num_c_dag = -1, num_c = -1;
    num_c_dag = rng(det_size);
    if (pauli_prob == 0.0) num_c = rng(det_size);

#ifdef EXT_DEBUG
    std::cerr << "* Proposing to remove: ";
    std::cerr << num_c_dag << "-th Cdag(" << block_index << ",...), ";
    std::cerr << num_c << "-th C(" << block_index << ",...)" << std::endl;
#endif

    // now mark 2 nodes for deletion
    if (pauli_prob == 0.0) tau1 = data.imp_trace.try_delete(num_c, block_index, false);
    tau2 = data.imp_trace.try_delete(num_c_dag, block_index, true);

    mc_weight_t t_ratio = 1.;

    if (pauli_prob > 0.0) {

      t_ratio = block_size * config.beta() / double(det_size);

      int rs_dag = det.get_x(num_c_dag).second;

      if (det_size > Nmax) {
        while (det_size > Nmax) Nmax *= 2;
        vec_ind.reserve(Nmax);
      }
      vec_ind.clear();

      // Look at position of tau2 among annihilation operators of the same flavor
      int j = 0, op_pos = -1;
      for (int i = 0; i < det_size; ++i) {
        if (det.get_y(i).first < tau2 && op_pos == -1) op_pos = j;
        if (det.get_y(i).second != rs_dag) continue;
        vec_ind.push_back(i);
        ++j;
      }

      int size = vec_ind.size();

      if (op_pos == -1) op_pos = size;

      int ic_nodagR = -1, ic_nodagL = -1, num_pauli = 0;

      // Check annihilation operators before and after tau2 (of the same flavor)
      if (size != 0) {
        ic_nodagR = (op_pos == size ? vec_ind[0] : vec_ind[op_pos]);
        ic_nodagL = (op_pos == 0 ? vec_ind[size-1] : vec_ind[op_pos-1]);
        num_pauli = (ic_nodagR == ic_nodagL ? 1 : 2);
      }

      if (num_pauli > 0) {
        bool pauli_move = (num_pauli == det_size);
        if (num_pauli != det_size) {
          double ran = rng();
          pauli_move = (ran <= pauli_prob);
        }
        if (pauli_move) { // Choose num_c between ic_nodagR and ic_nodagL
          if (num_pauli == 1)
            num_c = ic_nodagR;
          else {
            int ran_pauli = rng(2);
            num_c = (ran_pauli == 0 ? ic_nodagR : ic_nodagL);
          }
          if (num_pauli == det_size)
            t_ratio /= double(det_size);
          else
            t_ratio *= pauli_prob / double(num_pauli) + (1. - pauli_prob) / double(det_size);
        }
        else {  // Choose num_c uniformly
          num_c = rng(det_size);
          if (num_c == ic_nodagR || num_c == ic_nodagL)
            t_ratio *= pauli_prob / double(num_pauli) + (1. - pauli_prob) / double(det_size);
          else
            t_ratio *= (1. - pauli_prob) / double(det_size);
        }
      }
      else {
        num_c = rng(det_size);
        t_ratio /= double(det_size);
      }

      tau1 = data.imp_trace.try_delete(num_c, block_index, false);

      int rs = det.get_y(num_c).second;

      if (rs != rs_dag)
        t_ratio *= block_size * config.beta() / (1. - pauli_prob);
      else {
        if (size == 1)
          t_ratio /= (pauli_prob + (1. - pauli_prob) / double(block_size)) / config.beta();
        else {
          // Find closest annihilation operators to tau2 of the same flavor (excluding num_c)
          for (j = 0 ; j < size; ++j) {
            if (det.get_y(vec_ind[j]).first < tau2) break;
          }
          if (j == size) j = 0;
          if (vec_ind[j] == num_c) j = (j == size - 1 ? 0 : j + 1);
          ic_nodagR = vec_ind[j];
          j = (j == 0 ? size - 1 : j - 1);
          if (vec_ind[j] == num_c) j = (j == 0 ? size - 1 : j - 1);
          ic_nodagL = vec_ind[j];

          auto tR_nodag = det.get_y(ic_nodagR).first;
          auto tL_nodag = det.get_y(ic_nodagL).first;

          vec_ind.clear();

          // Find closest creation operators to tau2 of the same flavor
          j = 0; op_pos = -1;
          for (int i = 0; i < det_size; ++i) {
            if (det.get_x(i).first < tau2 && op_pos == -1) op_pos = j;
            if (det.get_x(i).second != rs_dag) continue;
            vec_ind.push_back(i);
            ++j;
          }

          size = vec_ind.size();

          if (size == 1)
            t_ratio /= (pauli_prob + (1. - pauli_prob) / double(block_size)) / config.beta();
          else {
            if (op_pos == -1) op_pos = size;
            int ic_dagR = (op_pos == size ? vec_ind[0] : vec_ind[op_pos]);
            op_pos = (op_pos == 0 ? size - 1 : op_pos - 1);
            op_pos = (op_pos == 0 ? size - 1 : op_pos - 1); // Go back two positions to skip num_c_dag
            int ic_dagL = vec_ind[op_pos];

            auto tR_dag = det.get_x(ic_dagR).first;
            auto tL_dag = det.get_x(ic_dagL).first;

            auto tR = ((tau2 - tR_dag) > (tau2 - tR_nodag) ? tR_nodag : tR_dag);
            auto tL = ((tL_dag - tau2) > (tL_nodag - tau2) ? tL_nodag : tL_dag);

            if (tR == tR_dag) {
              if ((tau2 - tau1) < (tau2 - tR))
                t_ratio /= pauli_prob / double(tau2 - tR) + (1. - pauli_prob) / (block_size * config.beta());
              else
                t_ratio *= block_size * config.beta() / (1. - pauli_prob);
            }
            else {
              if ((tau1 - tau2) < (tL - tau2))
                t_ratio /= pauli_prob / double(tL - tau2) + (1. - pauli_prob) / (block_size * config.beta());
              else
                t_ratio *= block_size * config.beta() / (1. - pauli_prob);
            }
          }
        }
      }
    }

    // record the length of the proposed removal
    dtau = double(tau2 - tau1);
    if (histo_proposed) *histo_proposed << dtau;

    auto det_ratio = det.try_remove(num_c_dag, num_c);

    // proposition probability
    if (pauli_prob == 0.0) t_ratio = std::pow(block_size * config.beta() / double(det_size), 2); // Size of the det before the try_delete!

    // For quick abandon
    double random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double p_yee = std::abs(det_ratio / t_ratio / data.atomic_weight);

    // recompute the atomic_weight
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) {
#ifdef EXT_DEBUG
      std::cerr << "atomic_weight == 0" << std::endl;
#endif
      return 0;
    }
    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(remove) atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                          << new_atomic_weight / data.atomic_weight << " in config " << config.get_id();

    mc_weight_t p = atomic_weight_ratio * det_ratio;

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

    if (!isfinite(p / t_ratio)){
      TRIQS_RUNTIME_ERROR << "(remove) p / t_ratio not finite p : " << p << " t_ratio :  " << t_ratio << " in config " << config.get_id();
    }
    return p / t_ratio;
  }

  mc_weight_t move_remove_c_cdag::accept() {

    // remove from the tree
    data.imp_trace.confirm_delete();

    // remove from the configuration
    config.erase(tau1);
    config.erase(tau2);
    config.finalize();

    // remove from the determinants
    data.dets[block_index].complete_operation();
    data.update_sign();
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;
    if (histo_accepted) *histo_accepted << dtau;

#ifdef EXT_DEBUG
    std::cerr << "* Move move_remove_c_cdag accepted" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    check_det_sequence(data.dets[block_index], config.get_id());
#endif

    return data.current_sign / data.old_sign;
  }

  void move_remove_c_cdag::reject() {

    config.finalize();
    data.imp_trace.cancel_delete();
    data.dets[block_index].reject_last_try();

#ifdef EXT_DEBUG
    std::cerr << "* Move move_remove_c_cdag rejected" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    check_det_sequence(data.dets[block_index], config.get_id());
#endif
  }
}
