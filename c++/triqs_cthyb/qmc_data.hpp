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
#include "impurity_trace.hpp"
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/det_manip.hpp>

namespace triqs_cthyb {
  using namespace triqs::gfs;
  using namespace triqs::mesh;
  using namespace nda;

  /************************
 * The Monte Carlo data
 ***********************/
  struct qmc_data {

    configuration config; // Configuration
    time_segment tau_seg;
    std::map<std::pair<int, int>, int> linindex; // Linear index constructed from block and inner indices
    atom_diag const &h_diag;                     // Diagonalization of the atomic problem
    mutable impurity_trace imp_trace;            // Calculator of the trace
    std::vector<int> n_inner;
    block_gf<imtime, delta_target_t> delta; // Hybridization function

    /// This callable object adapts the Delta function for the call of the det.
    struct delta_block_adaptor {
      gf<imtime, delta_target_t> delta_block; // make a copy. Needed in the real case anyway.

      delta_block_adaptor(gf_const_view<imtime, delta_target_t> delta_block) : delta_block(std::move(delta_block)) {}
      delta_block_adaptor(delta_block_adaptor const &) = default;
      delta_block_adaptor(delta_block_adaptor &&)      = default;
      delta_block_adaptor &operator=(delta_block_adaptor const &) = delete;
      delta_block_adaptor &operator=(delta_block_adaptor &&) = default;

      det_scalar_t operator()(std::pair<time_pt, int> const &x, std::pair<time_pt, int> const &y) const {
        det_scalar_t res = delta_block[closest_mesh_pt(double(x.first - y.first))](x.second, y.second);
        return (x.first >= y.first ? res : -res); // x,y first are time_pt, wrapping is automatic in the - operation, but need to
                                                  // compute the sign
      }

      friend void swap(delta_block_adaptor &dba1, delta_block_adaptor &dba2) noexcept { std::swap(dba1.delta_block, dba2.delta_block); }
    };

    std::vector<det_manip::det_manip<delta_block_adaptor>> dets; // The determinants
    int current_sign, old_sign;                                  // Permutation prefactor
    h_scalar_t atomic_weight;                                    // The current value of the trace or norm
    h_scalar_t atomic_reweighting;                               // The current value of the reweighting

    // Construction
    qmc_data(double beta, solve_parameters_t const &p, atom_diag const &h_diag, std::map<std::pair<int, int>, int> linindex,
             block_gf_const_view<imtime> delta, std::vector<int> n_inner, histo_map_t *histo_map)
       : config(beta),
         tau_seg(beta),
         linindex(linindex),
         h_diag(h_diag),
         imp_trace(beta, h_diag, histo_map, p.use_norm_as_weight, p.measure_density_matrix, p.performance_analysis),
         n_inner(n_inner),
         delta(map([](gf_const_view<imtime> d) { return real(d); }, delta)),
         current_sign(1),
         old_sign(1) {
      std::tie(atomic_weight, atomic_reweighting) = imp_trace.compute();
      dets.clear();
      for (auto const &bl : range(delta.size())) {
#ifdef HYBRIDISATION_IS_COMPLEX
        dets.emplace_back(delta_block_adaptor(delta[bl]), p.det_init_size);
#else
        if (!is_gf_real(delta[bl], 1e-10)) {
          //TRIQS_RUNTIME_ERROR << "The Delta(tau) block number " << bl << " is not real in tau space";
          if (p.verbosity >= 2) {
            std::cerr << "WARNING: The Delta(tau) block number " << bl << " is not real in tau space\n";
            std::cerr << "WARNING: max(Im[Delta(tau)]) = " << max_element(abs(imag(delta[bl].data()))) << "\n";
            std::cerr << "WARNING: Dissregarding the imaginary component in the calculation.\n";
          }
        }
        dets.emplace_back(delta_block_adaptor(real(delta[bl])), p.det_init_size);
#endif
        dets.back().set_singular_threshold(p.det_singular_threshold);
        dets.back().set_n_operations_before_check(p.det_n_operations_before_check);
        dets.back().set_precision_warning(p.det_precision_warning);
        dets.back().set_precision_error(p.det_precision_error);
      }
    }

    qmc_data(qmc_data const &) = delete; // Member imp_trace is not copyable
    qmc_data &operator=(qmc_data const &) = delete;

    void update_sign() {

      int s             = 0;
      size_t num_blocks = dets.size();
      std::vector<int> n_op_with_a_equal_to(num_blocks, 0), n_ndag_op_with_a_equal_to(num_blocks, 0);

      // In this first part we compute the sign to bring the configuration to
      // d^_1 d^_1 d^_1 ... d_1 d_1 d_1   d^_2 d^_2 ... d_2 d_2   ...   d^_n .. d_n

      // loop over the operators "op" in the trace (right to left)
      for (auto const &op : config) {

        // how many operators with an 'a' larger than "op" are there on the left of "op"?
        for (int a = op.second.block_index + 1; a < num_blocks; ++a) s += n_op_with_a_equal_to[a];
        n_op_with_a_equal_to[op.second.block_index]++;

        // if "op" is not a dagger how many operators of the same a but with a dagger are there on his right?
        if (op.second.dagger)
          s += n_ndag_op_with_a_equal_to[op.second.block_index];
        else
          n_ndag_op_with_a_equal_to[op.second.block_index]++;
      }

      // Now we compute the sign to bring the configuration to
      // d_1 d^_1 d_1 d^_1 ... d_1 d^_1   ...   d_n d^_n ... d_n d^_n
      for (int block_index = 0; block_index < num_blocks; block_index++) {
        int n = dets[block_index].size();
        s += n * (n + 1) / 2;
      }

      old_sign     = current_sign;
      current_sign = (s % 2 == 0 ? 1 : -1);
    }
  };

  //--------- DEBUG ---------

  using det_type = det_manip::det_manip<qmc_data::delta_block_adaptor>;

  // Print taus of operator sequence in dets
  inline void print_det_sequence(qmc_data const &data) {
    int i;
    int block_index;
    for (block_index = 0; block_index < data.dets.size(); ++block_index) {
      auto det = data.dets[block_index];
      if (det.size() == 0) return;
      std::cout << "BLOCK = " << block_index << std::endl;
      for (i = 0; i < det.size(); ++i) { // c_dag
        std::cout << " ic_dag = " << i << ": tau = " << det.get_x(i).first << std::endl;
      }
      for (i = 0; i < det.size(); ++i) { // c
        std::cout << " ic     = " << i << ": tau = " << det.get_y(i).first << std::endl;
      }
    }
  }

  inline void print_det_sequence(det_type const &det) {
    int i;
    for (i = 0; i < det.size(); ++i) { // c_dag
      std::cout << " ic_dag = " << i << ": tau = " << det.get_x(i).first << std::endl;
    }
    for (i = 0; i < det.size(); ++i) { // c
      std::cout << " ic     = " << i << ": tau = " << det.get_y(i).first << std::endl;
    }
  }

  // Check if dets are correctly ordered, otherwise complain
  inline void check_det_sequence(det_type const &det, int config_id) {
    if (det.size() == 0) return;
    auto tau = det.get_x(0).first;
    for (auto ic_dag = 0; ic_dag < det.size(); ++ic_dag) { // c_dag
      if (tau < det.get_x(ic_dag).first) {
        std::cout << "ERROR in det order in config " << config_id << std::endl;
        print_det_sequence(det);
        TRIQS_RUNTIME_ERROR << "Det ordering wrong: tau(ic_dag = " << ic_dag << ") = " << double(det.get_x(ic_dag).first);
      }
      tau = det.get_x(ic_dag).first;
    }
    tau = det.get_y(0).first;
    for (auto ic = 0; ic < det.size(); ++ic) { // c
      if (tau < det.get_y(ic).first) {
        std::cout << "ERROR in det order in config " << config_id << std::endl;
        print_det_sequence(det);
        TRIQS_RUNTIME_ERROR << "Det ordering wrong: tau(ic     = " << ic << ") = " << double(det.get_y(ic).first);
      }
      tau = det.get_y(ic).first;
    }
  }
} // namespace triqs_cthyb
