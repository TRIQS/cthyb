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
#include "./solver_core.hpp"
#include "./qmc_data.hpp"
#include <triqs/utility/callbacks.hpp>
#include <triqs/utility/exceptions.hpp>
#include <triqs/utility/variant_int_string.hpp>
#include <triqs/gfs.hpp>
#include <fstream>

#include "./moves/insert.hpp"
#include "./moves/remove.hpp"
#include "./moves/double_insert.hpp"
#include "./moves/double_remove.hpp"
#include "./moves/shift.hpp"
#include "./moves/global.hpp"
#include "./measures/g.hpp"
#include "./measures/g_legendre.hpp"
#include "./measures/perturbation_hist.hpp"
#include "./measures/density_matrix.hpp"
#include "./measures/average_sign.hpp"
#ifdef MEASURE_G2
#include "measure_g2.hpp"
#endif

namespace cthyb {

  struct index_visitor {
    std::vector<std::string> indices;
    void operator()(int i) { indices.push_back(std::to_string(i)); }
    void operator()(std::string s) { indices.push_back(s); }
  };

  solver_core::solver_core(double beta_, std::map<std::string, indices_type> const &gf_struct_, int n_iw, int n_tau, int n_l)
     : beta(beta_), gf_struct(gf_struct_) {

    if (n_tau < 2 * n_iw)
      TRIQS_RUNTIME_ERROR << "Must use as least twice as many tau points as Matsubara frequencies: n_iw = " << n_iw << " but n_tau = " << n_tau
                          << ".";

    std::vector<std::string> block_names;
    std::vector<gf<imfreq>> g0_iw_blocks;
    std::vector<gf<imtime>> g_tau_blocks;
    std::vector<gf<legendre>> g_l_blocks;
    std::vector<gf<imtime>> delta_tau_blocks;
    std::vector<gf<imtime, delta_target_t>> g_tau_accum_blocks; //  Local real or complex (if complex mode) quantities for accumulation

    for (auto const &bl : gf_struct) {
      block_names.push_back(bl.first);
      int n = bl.second.size();

      index_visitor iv;
      for (auto &ind : bl.second) { apply_visitor(iv, ind); }
      std::vector<std::vector<std::string>> indices{{iv.indices, iv.indices}};

      g0_iw_blocks.push_back(gf<imfreq>{{beta, Fermion, n_iw}, {n, n}, indices});
      g_tau_blocks.push_back(gf<imtime>{{beta, Fermion, n_tau}, {n, n}, indices});
      g_l_blocks.push_back(gf<legendre>{{beta, Fermion, static_cast<size_t>(n_l)}, {n, n}, indices}); // FIXME: cast is ugly
      delta_tau_blocks.push_back(gf<imtime>{{beta, Fermion, n_tau}, {n, n}, indices});
      g_tau_accum_blocks.push_back(gf<imtime, delta_target_t>{{beta, Fermion, n_tau}, {n, n}});
    }

    _G0_iw       = make_block_gf(block_names, g0_iw_blocks);
    _G_tau       = make_block_gf(block_names, g_tau_blocks);
    _G_l         = make_block_gf(block_names, g_l_blocks);
    _Delta_tau   = make_block_gf(block_names, delta_tau_blocks);
    _G_tau_accum = make_block_gf(block_names, g_tau_accum_blocks);

#ifdef MEASURE_G2
    // Fill G^2 containers with empty blocks
    using g2_inu_v_t = std::vector<g2_iw_inu_inup_block_t>;
    std::vector<g2_inu_v_t> g2_inu_blocks(gf_struct.size(), g2_inu_v_t(gf_struct.size()));
    _G2_iw_inu_inup_pp    = make_block2_gf(block_names, block_names, g2_inu_blocks);
    _G2_iw_inu_inup_ph    = make_block2_gf(block_names, block_names, g2_inu_blocks);
    using g2_legendre_v_t = std::vector<g2_iw_l_lp_block_t>;
    std::vector<g2_legendre_v_t> g2_legendre_blocks(gf_struct.size(), g2_legendre_v_t(gf_struct.size()));
    _G2_iw_l_lp_pp = make_block2_gf(block_names, block_names, g2_legendre_blocks);
    _G2_iw_l_lp_ph = make_block2_gf(block_names, block_names, g2_legendre_blocks);
#endif
  }

  /// -------------------------------------------------------------------------------------------

  void solver_core::solve(solve_parameters_t const &params) {

    _last_solve_parameters = params;

    // determine basis of operators to use
    fundamental_operator_set fops;
    for (auto const &bl : gf_struct) {
      for (auto const &a : bl.second) { fops.insert(bl.first, a); }
    }

    // setup the linear index map
    std::map<std::pair<int, int>, int> linindex;
    int block_index = 0;
    for (auto const &bl : gf_struct) {
      int inner_index = 0;
      for (auto const &a : bl.second) {
        linindex[std::make_pair(block_index, inner_index)] = fops[{bl.first, a}];
        inner_index++;
      }
      block_index++;
    }

    // Make list of block sizes
    std::vector<int> n_inner;
    for (auto const &bl : gf_struct) { n_inner.push_back(bl.second.size()); }

    // Calculate imfreq quantities
    auto G0_iw_inv = map([](gf_const_view<imfreq> x) { return triqs::gfs::inverse(x); }, _G0_iw);
    auto Delta_iw  = G0_iw_inv;

    _h_loc = params.h_int;

    // Do I have imaginary components in my local Hamiltonian?
    auto max_imag                                  = 0.0;
    for (int b : range(gf_struct.size())) max_imag = std::max(max_imag, max_element(abs(real(_G0_iw[b].singularity()(2)))));

    // Add quadratic terms to h_loc
    int b = 0;
    for (auto const &bl : gf_struct) {
      int n1 = 0;
      for (auto const &a1 : bl.second) {
        int n2 = 0;
        for (auto const &a2 : bl.second) {
#ifdef LOCAL_HAMILTONIAN_IS_COMPLEX
          dcomplex e_ij;
          if (max_imag > params.imag_threshold)
            e_ij = _G0_iw[b].singularity()(2)(n1, n2);
          else
            e_ij = _G0_iw[b].singularity()(2)(n1, n2).real();
#else
          auto e_ij = _G0_iw[b].singularity()(2)(n1, n2).real();
#endif
          _h_loc = _h_loc + e_ij * c_dag<h_scalar_t>(bl.first, a1) * c<h_scalar_t>(bl.first, a2);
          n2++;
        }
        n1++;
      }
      b++;
    }

    // Determine terms Delta_iw from G0_iw and ensure that the 1/iw behaviour of G0_iw is correct
    b = 0;
    range _;
    triqs::clef::placeholder<0> iw_;
    for (auto const &bl : gf_struct) {
      Delta_iw[b](iw_) << G0_iw_inv[b].singularity()(-1) * iw_ + G0_iw_inv[b].singularity()(0);
      Delta_iw[b]     = Delta_iw[b] - G0_iw_inv[b];
      _Delta_tau[b]() = inverse_fourier(Delta_iw[b]);
      // Force all diagonal elements to be real
      for (int i : range(bl.second.size())) _Delta_tau[b].data()(_, i, i) = real(_Delta_tau[b].data()(_, i, i));
      // If off-diagonal elements are below threshold, set to real
      if (max_element(abs(imag(_Delta_tau[b].data()))) < params.imag_threshold) _Delta_tau[b].data() = real(_Delta_tau[b].data());
      b++;
    }

    // Report what h_loc we are using
    if (params.verbosity >= 2) std::cout << "The local Hamiltonian of the problem:" << std::endl << _h_loc << std::endl;
    // Reset the histograms
    _performance_analysis.clear();
    histo_map_t *histo_map = params.performance_analysis ? &_performance_analysis : nullptr;

    // Determine block structure
    if (params.partition_method == "autopartition") {
      if (params.verbosity >= 2) std::cout << "Using autopartition algorithm to partition the local Hilbert space" << std::endl;
      h_diag = {_h_loc, fops};
    } else if (params.partition_method == "quantum_numbers") {
      if (params.quantum_numbers.empty()) TRIQS_RUNTIME_ERROR << "No quantum numbers provided.";
      if (params.verbosity >= 2) std::cout << "Using quantum numbers to partition the local Hilbert space" << std::endl;
      h_diag = {_h_loc, fops, params.quantum_numbers};
    } else if (params.partition_method == "none") { // give empty quantum numbers list
      std::cout << "Not partitioning the local Hilbert space" << std::endl;
      h_diag = {_h_loc, fops, std::vector<many_body_op_t>{}};
    } else
      TRIQS_RUNTIME_ERROR << "Partition method " << params.partition_method << " not recognised.";

    // FIXME save h_loc to be able to rebuild h_diag in an analysis program.
    //if (_comm.rank() ==0) h5_write(h5::file("h_loc.h5",'w'), "h_loc", _h_loc, fops);

    if (params.verbosity >= 2) std::cout << "Found " << h_diag.n_blocks() << " subspaces." << std::endl;

    if (params.performance_analysis) std::ofstream("impurity_blocks.dat") << h_diag;

    // If one is interested only in the atomic problem
    if (params.n_warmup_cycles == 0 && params.n_cycles == 0) {
      if (params.measure_density_matrix) _density_matrix = atomic_density_matrix(h_diag, beta);
      return;
    }

    // Initialise Monte Carlo quantities
    qmc_data data(beta, params, h_diag, linindex, _Delta_tau, n_inner, histo_map);
    auto qmc = mc_tools::mc_generic<mc_weight_t>(params.random_name, params.random_seed, 1.0, params.verbosity);

    // Moves
    using move_set_type = mc_tools::move_set<mc_weight_t>;
    move_set_type inserts(qmc.get_rng());
    move_set_type removes(qmc.get_rng());
    move_set_type double_inserts(qmc.get_rng());
    move_set_type double_removes(qmc.get_rng());

    auto &delta_names  = _Delta_tau.block_names();
    auto get_prob_prop = [&params](std::string const &block_name) {
      auto f = params.proposal_prob.find(block_name);
      return (f != params.proposal_prob.end() ? f->second : 1.0);
    };
    for (size_t block = 0; block < _Delta_tau.size(); ++block) {
      int block_size         = _Delta_tau[block].data().shape()[1];
      auto const &block_name = delta_names[block];
      double prop_prob       = get_prob_prop(block_name);
      inserts.add(move_insert_c_cdag(block, block_size, block_name, data, qmc.get_rng(), histo_map), "Insert Delta_" + block_name, prop_prob);
      removes.add(move_remove_c_cdag(block, block_size, block_name, data, qmc.get_rng(), histo_map), "Remove Delta_" + block_name, prop_prob);
      if (params.move_double) {
        for (size_t block2 = 0; block2 < _Delta_tau.size(); ++block2) {
          int block_size2         = _Delta_tau[block2].data().shape()[1];
          auto const &block_name2 = delta_names[block2];
          double prop_prob2       = get_prob_prop(block_name2);
          double_inserts.add(
             move_insert_c_c_cdag_cdag(block, block2, block_size, block_size2, block_name, block_name2, data, qmc.get_rng(), histo_map),
             "Insert Delta_" + block_name + "_" + block_name2, prop_prob * prop_prob2);
          double_removes.add(
             move_remove_c_c_cdag_cdag(block, block2, block_size, block_size2, block_name, block_name2, data, qmc.get_rng(), histo_map),
             "Remove Delta_" + block_name + "_" + block_name2, prop_prob * prop_prob2);
        }
      }
    }

    qmc.add_move(std::move(inserts), "Insert two operators", 1.0);
    qmc.add_move(std::move(removes), "Remove two operators", 1.0);
    if (params.move_double) {
      qmc.add_move(std::move(double_inserts), "Insert four operators", 1.0);
      qmc.add_move(std::move(double_removes), "Remove four operators", 1.0);
    }

    if (params.move_shift) qmc.add_move(move_shift_operator(data, qmc.get_rng(), histo_map), "Shift one operator", 1.0);

    if (params.move_global.size()) {
      move_set_type global(qmc.get_rng());
      for (auto const &mv : params.move_global) {
        auto const &name          = mv.first;
        auto const &substitutions = mv.second;
        global.add(move_global(name, substitutions, data, qmc.get_rng()), name, 1.0);
      }
      qmc.add_move(std::move(global), "Global moves", params.move_global_prob);
    }

#ifdef MEASURE_G2
    if (params.measure_g2_inu || params.measure_g2_legendre) {
      auto g2_blocks_to_measure = params.measure_g2_blocks;

      // Measure all blocks
      if (g2_blocks_to_measure.empty()) {
        for (auto const &bn1 : gf_struct) {
          for (auto const &bn2 : gf_struct) { g2_blocks_to_measure.emplace(bn1.first, bn2.first); }
        }
      } else { // Check the blocks we've been asked to measure
        for (auto const &bn : g2_blocks_to_measure) {
          if (!gf_struct.count(bn.first)) TRIQS_RUNTIME_ERROR << "Invalid left block name " << bn.first << " for G^2 measurement";
          if (!gf_struct.count(bn.second)) TRIQS_RUNTIME_ERROR << "Invalid right block name " << bn.second << " for G^2 measurement";
        }
      }

      if (!params.measure_g2_pp && !params.measure_g2_ph)
        TRIQS_RUNTIME_ERROR << "You must switch on at least one of measure_g2_pp and measure_g2_ph!";

      auto const &delta_names = _Delta_tau.domain().names();
      for (int b1 = 0; b1 < delta_names.size(); ++b1) {
        for (int b2 = 0; b2 < delta_names.size(); ++b2) {
          auto const &bn1 = delta_names[b1];
          auto const &bn2 = delta_names[b2];
          if (!g2_blocks_to_measure.count({bn1, bn2})) continue;

          int s1 = get_target_shape(_Delta_tau[b1])[0];
          int s3 = get_target_shape(_Delta_tau[b2])[0];
          int s2 = params.measure_g2_block_order == AABB ? s1 : s3;
          int s4 = params.measure_g2_block_order == AABB ? s3 : s1;

          auto make_measure_name = [&bn1, &bn2](bool legendre, g2_channel channel, block_order bo) {
            return std::string("G^2 measure, ") + (legendre ? "Legendre" : "Matsubara") + ", " + (channel == PP ? "pp" : "ph")
               + (bo == AABB ? (" (" + bn1 + "," + bn1 + "," + bn2 + "," + bn2 + ")") : (" (" + bn1 + "," + bn2 + "," + bn2 + "," + bn1 + ")"));
          };

          int n_iw = params.measure_g2_n_iw;

          // Matsubara measurements
          if (params.measure_g2_inu) {
            int n_inu = params.measure_g2_n_inu;
            if (params.measure_g2_pp) {
              auto &block = _G2_iw_inu_inup_pp[{b1, b2}];
              block       = g2_iw_inu_inup_block_t{{{beta, Boson, n_iw}, {beta, Fermion, n_inu}, {beta, Fermion, n_inu}}, {s1, s2, s3, s4}};
              if (params.measure_g2_block_order == AABB)
                qmc.add_measure(measure_g2_inu<PP, AABB>(b1, b2, block, data), make_measure_name(false, PP, AABB));
              else
                qmc.add_measure(measure_g2_inu<PP, ABBA>(b1, b2, block, data), make_measure_name(false, PP, ABBA));
            }

            if (params.measure_g2_ph) {
              auto &block = _G2_iw_inu_inup_ph[{b1, b2}];
              block       = g2_iw_inu_inup_block_t{{{beta, Boson, n_iw}, {beta, Fermion, n_inu}, {beta, Fermion, n_inu}}, {s1, s2, s3, s4}};
              if (params.measure_g2_block_order == AABB)
                qmc.add_measure(measure_g2_inu<PH, AABB>(b1, b2, block, data), make_measure_name(false, PH, AABB));
              else
                qmc.add_measure(measure_g2_inu<PH, ABBA>(b1, b2, block, data), make_measure_name(false, PH, ABBA));
            }
          }

          // Legendre measurements
          if (params.measure_g2_legendre) {
            size_t n_l = params.measure_g2_n_l;
            if (params.measure_g2_pp) {
              auto &block = _G2_iw_l_lp_pp[{b1, b2}];
              block       = g2_iw_l_lp_block_t{{{beta, Boson, n_iw}, {beta, Fermion, n_l}, {beta, Fermion, n_l}}, {s1, s2, s3, s4}};
              if (params.measure_g2_block_order == AABB)
                qmc.add_measure(measure_g2_legendre<PP, AABB>(b1, b2, block, data), make_measure_name(true, PP, AABB));
              else
                qmc.add_measure(measure_g2_legendre<PP, ABBA>(b1, b2, block, data), make_measure_name(true, PP, ABBA));
            }

            if (params.measure_g2_ph) {
              auto &block = _G2_iw_l_lp_ph[{b1, b2}];
              block       = g2_iw_l_lp_block_t{{{beta, Boson, n_iw}, {beta, Fermion, n_l}, {beta, Fermion, n_l}}, {s1, s2, s3, s4}};
              if (params.measure_g2_block_order == AABB)
                qmc.add_measure(measure_g2_legendre<PH, AABB>(b1, b2, block, data), make_measure_name(true, PH, AABB));
              else
                qmc.add_measure(measure_g2_legendre<PH, ABBA>(b1, b2, block, data), make_measure_name(true, PH, ABBA));
            }
          }
        }
      }
    }
#endif

    // Measurements
    if (params.measure_g_tau) {
      auto &g_names = _G_tau.block_names();
      for (size_t block = 0; block < _G_tau.size(); ++block) {
        qmc.add_measure(measure_g(block, _G_tau_accum[block], data), "G measure (" + g_names[block] + ")");
      }
    }
    if (params.measure_g_l) {
      auto &g_names = _G_l.block_names();
      for (size_t block = 0; block < _G_l.size(); ++block) {
        qmc.add_measure(measure_g_legendre(block, _G_l[block], data), "G_l measure (" + g_names[block] + ")");
      }
    }
    if (params.measure_pert_order) {
      auto &g_names = _G_tau.block_names();
      for (size_t block = 0; block < _G_tau.size(); ++block) {
        auto const &block_name = g_names[block];
        qmc.add_measure(measure_perturbation_hist(block, data, _pert_order[block_name]), "Perturbation order (" + block_name + ")");
      }
      qmc.add_measure(measure_perturbation_hist_total(data, _pert_order_total), "Perturbation order");
    }
    if (params.measure_density_matrix) {
      if (!params.use_norm_as_weight)
        TRIQS_RUNTIME_ERROR << "To measure the density_matrix of atomic states, you need to set "
                               "use_norm_as_weight to True, i.e. to reweight the QMC";
      qmc.add_measure(measure_density_matrix{data, _density_matrix}, "Density Matrix for local static observable");
    }

    qmc.add_measure(measure_average_sign{data, _average_sign}, "Average sign");

    // Run! The empty (starting) configuration has sign = 1
    _solve_status =
       qmc.warmup_and_accumulate(params.n_warmup_cycles, params.n_cycles, params.length_cycle, triqs::utility::clock_callback(params.max_time));
    qmc.collect_results(_comm);

    if (params.verbosity >= 2) std::cout << "Average sign: " << _average_sign << std::endl;

    // Copy local (real or complex) G_tau back to complex G_tau
    if (params.measure_g_tau) _G_tau = _G_tau_accum;
  }
}
