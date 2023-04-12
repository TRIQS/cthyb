/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
 * Copyright (C) 2017, H. UR Strand, P. Seth, I. Krivenko,
 *                     M. Ferrero and O. Parcollet
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
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <fstream>
#include <variant>

#include "./moves/insert.hpp"
#include "./moves/remove.hpp"
#include "./moves/double_insert.hpp"
#include "./moves/double_remove.hpp"
#include "./moves/shift.hpp"
#include "./moves/global.hpp"
#include "./measures/G_tau.hpp"
#include "./measures/G_l.hpp"
#include "./measures/O_tau_ins.hpp"
#include "./measures/perturbation_hist.hpp"
#include "./measures/density_matrix.hpp"
#include "./measures/average_sign.hpp"
#include "./measures/average_order.hpp"
#include "./measures/auto_corr_time.hpp"
#ifdef CTHYB_G2_NFFT
#include "./measures/G2_tau.hpp"
#include "./measures/G2_iw.hpp"
#include "./measures/G2_iw_nfft.hpp"
#include "./measures/G2_iwll.hpp"
#endif
#include "./measures/util.hpp"

namespace triqs_cthyb {

  struct index_visitor {
    std::vector<std::string> indices;
    void operator()(int i) { indices.push_back(std::to_string(i)); }
    void operator()(std::string s) { indices.push_back(s); }
  };

  solver_core::solver_core(constr_parameters_t const &p)
     : beta(p.beta), gf_struct(p.gf_struct), n_iw(p.n_iw), n_tau(p.n_tau), n_l(p.n_l), delta_interface(p.delta_interface), constr_parameters(p) {

    if (p.n_tau < 2 * p.n_iw)
      TRIQS_RUNTIME_ERROR
         << "Must use as least twice as many tau points as Matsubara frequencies: n_iw = " << p.n_iw
         << " but n_tau = " << p.n_tau << ".";

    // Allocate single particle greens functions
    if (not delta_interface) _G0_iw = block_gf<imfreq>({beta, Fermion, n_iw}, gf_struct);
    _Delta_tau = block_gf<imtime>({beta, Fermion, n_tau}, gf_struct);
  }

  /// -------------------------------------------------------------------------------------------

  void solver_core::solve(solve_parameters_t const &solve_parameters_) {

    solve_parameters = solve_parameters_;
    solve_parameters_t params(solve_parameters_);

    // Merge constr_params and solve_params
    //params_t params(constr_parameters, solve_parameters);

    // http://patorjk.com/software/taag/#p=display&f=Calvin%20S&t=TRIQS%20cthyb
    if (params.verbosity >= 2)
      std::cout << "\n"
                   "╔╦╗╦═╗╦╔═╗ ╔═╗  ┌─┐┌┬┐┬ ┬┬ ┬┌┐ \n"
                   " ║ ╠╦╝║║═╬╗╚═╗  │   │ ├─┤└┬┘├┴┐\n"
                   " ╩ ╩╚═╩╚═╝╚╚═╝  └─┘ ┴ ┴ ┴ ┴ └─┘\n\n";

    // determine basis of operators to use
    fundamental_operator_set fops;
    for (auto const &[bl, bl_size] : gf_struct) {
      for (auto idx : range(bl_size)) { fops.insert(bl, idx); }
    }

    // setup the linear index map
    std::map<std::pair<int, int>, int> linindex;
    int block_index = 0;
    for (auto const &[bl, bl_size] : gf_struct) {
      int inner_index = 0;
      for (auto idx : range(bl_size)) {
        linindex[std::make_pair(block_index, inner_index)] = fops[{bl, idx}];
        inner_index++;
      }
      block_index++;
    }

    // Make list of block sizes
    std::vector<int> n_inner;
    for (auto const &[bl, bl_size] : gf_struct) { n_inner.push_back(bl_size); }

    if (not delta_interface) {

      // ==== Assert that G0_iw fulfills the fundamental property G(iw)[i,j] = G(-iw)*[j,i] ====

      if (not is_gf_hermitian(_G0_iw.value())) {
        if (params.verbosity >= 2)
          std::cout << "!-------------------------------------------------------------------------------------------!\n"
                       "! WARNING: S.G0_iw violates fundamental Green Function property G0(iw)[i,j] = G0(-iw)*[j,i] !\n"
                       "! Symmetrizing S.G0_iw ...                                                                  !\n"
                       "!-------------------------------------------------------------------------------------------!\n\n";
        _G0_iw = make_hermitian(_G0_iw.value());
      }

      // ==== Compute Delta from G0_iw ====

      // Initialize Delta with iw - inv(G0[iw])
      auto Delta_iw = inverse(_G0_iw.value());
      for (auto bl : range(gf_struct.size()))
        for (auto iw : Delta_iw[bl].mesh()) Delta_iw[bl][iw] = iw - Delta_iw[bl][iw];

      // Compute the constant part of Delta
      Delta_infty_vec = map(
         // Compute 0th moment of one block
         [imag_threshold = params.imag_threshold](gf_const_view<imfreq> d) {
           auto [tail, err] = fit_hermitian_tail(d);
           if (err > 1e-5) std::cerr << "WARNING: Tail fit in extraction of delta(infty) has large error of: " << err << std::endl;
           auto Delta_infty = matrix<dcomplex>{tail(0, ellipsis())};
#ifndef HYBRIDISATION_IS_COMPLEX
           double max_imag = max_element(abs(imag(Delta_infty)));
           if (max_imag > imag_threshold)
             TRIQS_RUNTIME_ERROR << "Largest imaginary element of delta(infty) e.g. of the local part of G0: " << max_imag
                                 << ", is larger than the set parameter imag_threshold " << imag_threshold;
#endif
           return Delta_infty;
         },
         Delta_iw);

      // Subtract constant part from Delta and perform Fourier transform
      for (auto bl : range(gf_struct.size())) {
        for (auto iw : Delta_iw[bl].mesh()) Delta_iw[bl][iw] = Delta_iw[bl][iw] - Delta_infty_vec.value()[bl];
        auto [Delta_tail_b, tail_err] = fit_hermitian_tail(Delta_iw[bl]);
        _Delta_tau[bl]()              = fourier(Delta_iw[bl], Delta_tail_b);
      }

      // ==== Compute h_loc ====

      _h_loc0 = {}; 

      // Add non-interacting terms to h_loc
      for (auto bl : range(gf_struct.size())) {
        for (auto [n1, n2] : Delta_iw[bl].target_indices()) {
#ifdef LOCAL_HAMILTONIAN_IS_COMPLEX
          dcomplex e_ij;
          double max_imag = max_element(abs(imag(Delta_infty_vec.value()[bl])));
          if (max_imag > params.imag_threshold)
            e_ij = Delta_infty_vec.value()[bl](n1, n2);
          else
            e_ij = Delta_infty_vec.value()[bl](n1, n2).real();
#else
          double e_ij = Delta_infty_vec.value()[bl](n1, n2).real();
#endif
          // Set off diagonal terms to 0 if they are below off_diag_threshold
          if (n1 != n2 && abs(Delta_infty_vec.value()[bl](n1, n2)) < params.off_diag_threshold) e_ij = 0.0;

	  auto bl_name = gf_struct[bl].first;
          _h_loc0 = _h_loc0 + e_ij * c_dag<h_scalar_t>(bl_name, n1) * c<h_scalar_t>(bl_name, n2);
        }
      }

    } else { // Delta Interface
      if (not is_gf_hermitian(_Delta_tau)) {
        if (params.verbosity >= 2)
          std::cout << "!---------------------------------------------------------------------------------------!\n"
                       "! WARNING: S.Delta_tau violates fundamental symmetry Delta(tau)[i,j] = Delta(tau)*[j,i] !\n"
                       "! Symmetrizing S.Delta_tau ...                                                          !\n"
                       "!---------------------------------------------------------------------------------------!\n\n";
        _Delta_tau = make_hermitian(_Delta_tau);
      }
      if (not params.h_loc0.has_value()) TRIQS_RUNTIME_ERROR << "h_loc0 must be provided when using the Delta interface";
      _h_loc0 = params.h_loc0.value();
    }

    _h_loc  = params.h_int + _h_loc0;

#ifndef HYBRIDISATION_IS_COMPLEX
    // Check that diagonal components of Delta_tau are real
    for (auto bl : range(gf_struct.size())) {
      long bl_size = gf_struct[bl].second;
      // Force all diagonal elements to be real
      auto _ = range::all;
      for (auto [i, j] : product_range(bl_size, bl_size)) {
        auto Delta_tau_bl_ij = _Delta_tau[bl].data()(_, i, j);
        double max_imag      = max_element(abs(imag(Delta_tau_bl_ij)));
        if (i == j && max_imag > 1e-10) {
          std::cout << "Warning! Delta_tau diagonal term has max imaginary part: " << max_imag << "Disregarding imaginary part \n";
          Delta_tau_bl_ij = real(Delta_tau_bl_ij);
        } else if (max_imag < params.imag_threshold) {
          Delta_tau_bl_ij = real(Delta_tau_bl_ij);
        }
      }
    }
#endif

    // Report what h_loc we are using
    if (params.verbosity >= 2)
      std::cout << "The local Hamiltonian of the problem:" << std::endl << _h_loc << std::endl;
    // Reset the histograms
    _performance_analysis.clear();
    histo_map_t *histo_map = params.performance_analysis ? &_performance_analysis : nullptr;

    // Determine block structure
    if (params.partition_method == "autopartition") {
      if (params.verbosity >= 2)
        std::cout << "Using autopartition algorithm to partition the local Hilbert space"
                  << std::endl;
      if (params.loc_n_min == 0 && params.loc_n_max == INT_MAX)
        h_diag = {_h_loc, fops};
      else {
        if (params.verbosity >= 2)
          std::cout << "Restricting the local Hilbert space to states with [" << params.loc_n_min
                    << ";" << params.loc_n_max << "] particles" << std::endl;
        h_diag = {_h_loc, fops, params.loc_n_min, params.loc_n_max};
      }
    } else if (params.partition_method == "quantum_numbers") {
      if (params.quantum_numbers.empty()) TRIQS_RUNTIME_ERROR << "No quantum numbers provided.";
      if (params.verbosity >= 2)
        std::cout << "Using quantum numbers to partition the local Hilbert space" << std::endl;
      h_diag = {_h_loc, fops, params.quantum_numbers};
    } else if (params.partition_method == "none") { // give empty quantum numbers list
      std::cout << "Not partitioning the local Hilbert space" << std::endl;
      h_diag = {_h_loc, fops, std::vector<many_body_op_t>{}};
    } else
      TRIQS_RUNTIME_ERROR << "Partition method " << params.partition_method << " not recognised.";

    // FIXME save h_loc to be able to rebuild h_diag in an analysis program.
    //if (_comm.rank() ==0) h5_write(h5::file("h_loc.h5",'w'), "h_loc", _h_loc, fops);

    if (params.verbosity >= 2)
      std::cout << "Found " << h_diag.n_subspaces() << " subspaces." << std::endl;

    if (params.performance_analysis) std::ofstream("impurity_blocks.dat") << h_diag;

    // If one is interested only in the atomic problem
    if (params.n_warmup_cycles == 0 && params.n_cycles == 0) {
      if (params.measure_density_matrix) _density_matrix = atomic_density_matrix(h_diag, beta);
      return;
    }

    // Initialise Monte Carlo quantities
    qmc_data data(beta, params, h_diag, linindex, _Delta_tau, n_inner, histo_map);
    auto qmc =
       mc_tools::mc_generic<mc_weight_t>(params.random_name, params.random_seed, params.verbosity);

    // --------------------------------------------------------------------------
    // Moves
    // --------------------------------------------------------------------------

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
      inserts.add(move_insert_c_cdag(block, block_size, block_name, data, qmc.get_rng(), histo_map),
                  "Insert Delta_" + block_name, prop_prob);
      removes.add(move_remove_c_cdag(block, block_size, block_name, data, qmc.get_rng(), histo_map),
                  "Remove Delta_" + block_name, prop_prob);
      if (params.move_double) {
        for (size_t block2 = 0; block2 < _Delta_tau.size(); ++block2) {
          int block_size2         = _Delta_tau[block2].data().shape()[1];
          auto const &block_name2 = delta_names[block2];
          double prop_prob2       = get_prob_prop(block_name2);
          double_inserts.add(
             move_insert_c_c_cdag_cdag(block, block2, block_size, block_size2, block_name,
                                       block_name2, data, qmc.get_rng(), histo_map),
             "Insert Delta_" + block_name + "_" + block_name2, prop_prob * prop_prob2);
          double_removes.add(
             move_remove_c_c_cdag_cdag(block, block2, block_size, block_size2, block_name,
                                       block_name2, data, qmc.get_rng(), histo_map),
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

    if (params.move_shift)
      qmc.add_move(move_shift_operator(data, qmc.get_rng(), histo_map), "Shift one operator", 1.0);

    if (params.move_global.size()) {
      move_set_type global(qmc.get_rng());
      for (auto const &mv : params.move_global) {
        auto const &name          = mv.first;
        auto const &substitutions = mv.second;
        global.add(move_global(name, substitutions, data, qmc.get_rng()), name, 1.0);
      }
      qmc.add_move(std::move(global), "Global moves", params.move_global_prob);
    }

    // --------------------------------------------------------------------------
    // Measurements
    // --------------------------------------------------------------------------

    // --------------------------------------------------------------------------
    // Two-particle correlators

    G2_measures_t G2_measures(_Delta_tau, gf_struct, params);

#ifdef CTHYB_G2_NFFT
    // Imaginary-time binning
    if (params.measure_G2_tau)
      qmc.add_measure(measure_G2_tau{G2_tau, data, G2_measures},
                      "G2_tau imaginary-time measurement");

    // NFFT Matsubara frequency measures

    if (params.measure_G2_iw_nfft)
      qmc.add_measure(measure_G2_iw_nfft<G2_channel::AllFermionic>{G2_iw_nfft, data, G2_measures},
                      "G2_iw nfft fermionic measurement");
    if (params.measure_G2_iw_pp_nfft)
      qmc.add_measure(measure_G2_iw_nfft<G2_channel::PP>{G2_iw_pp_nfft, data, G2_measures},
                      "G2_iw_pp nfft particle-particle measurement");
    if (params.measure_G2_iw_ph_nfft)
      qmc.add_measure(measure_G2_iw_nfft<G2_channel::PH>{G2_iw_ph_nfft, data, G2_measures},
                      "G2_iw_ph nfft particle-hole measurement");

    // Direct Matsubara frequency measurement

    if (params.measure_G2_iw)
      qmc.add_measure(measure_G2_iw<G2_channel::AllFermionic>{G2_iw, data, G2_measures},
                      "G2_iw fermionic measurement");
    if (params.measure_G2_iw_pp)
      qmc.add_measure(measure_G2_iw<G2_channel::PP>{G2_iw_pp, data, G2_measures},
                      "G2_iw_pp particle-particle measurement");
    if (params.measure_G2_iw_ph)
      qmc.add_measure(measure_G2_iw<G2_channel::PH>{G2_iw_ph, data, G2_measures},
                      "G2_iw_ph particle-hole measurement");

    // Legendre mixed basis measurements
    if (params.measure_G2_iwll_pp)
      qmc.add_measure(measure_G2_iwll<G2_channel::PP>{G2_iwll_pp, data, G2_measures},
                      "G2_iwll_pp Legendre particle-particle measurement");
    if (params.measure_G2_iwll_ph)
      qmc.add_measure(measure_G2_iwll<G2_channel::PH>{G2_iwll_ph, data, G2_measures},
                      "G2_iwll_ph Legendre particle-hole measurement");
#endif

    // --------------------------------------------------------------------------
    // Single-particle correlators

    if (params.measure_O_tau) {

      const auto &[O1, O2] = *params.measure_O_tau;
      auto comm_0          = O1 * O2 - O2 * O1;
      auto comm_1          = O1 * _h_loc - _h_loc * O1;
      auto comm_2          = O2 * _h_loc - _h_loc * O2;

      if (!comm_0.is_zero() || !comm_1.is_zero() || !comm_2.is_zero()) {
        if (params.verbosity >= 2) {
          TRIQS_RUNTIME_ERROR << "Error: measure_O_tau, supplied operators does not commute with "
                                 "the local Hamiltonian.\n"
                              << "[O1, O2] = " << comm_0 << "\n"
                              << "[O1, H_loc] = " << comm_1 << "\n"
                              << "[O2, H_loc] = " << comm_2 << "\n";
        }
      }
      qmc.add_measure(
         measure_O_tau_ins{O_tau, data, n_tau, O1, O2, params.measure_O_tau_min_ins, qmc.get_rng()},
         "O_tau insertion measure");
    }

    if (params.measure_G_tau) {
      G_tau = block_gf<imtime>{{beta, Fermion, n_tau}, gf_struct};
      qmc.add_measure(measure_G_tau{data, n_tau, gf_struct, container_set()}, "G_tau measure");
    }

    if (params.measure_G_l) qmc.add_measure(measure_G_l{G_l, data, n_l, gf_struct}, "G_l measure");

    // Other measurements
    if (params.measure_pert_order) {
      auto &g_names = _Delta_tau.block_names();
      perturbation_order       = histo_map_t{};
      perturbation_order_total = histogram{};
      for (size_t block = 0; block < _Delta_tau.size(); ++block) {
        auto const &block_name = g_names[block];
        qmc.add_measure(measure_perturbation_hist(block, data, (*perturbation_order)[block_name]), "Perturbation order (" + block_name + ")");
      }
      qmc.add_measure(measure_perturbation_hist_total(data, *perturbation_order_total), "Perturbation order");
    }
    if (params.measure_density_matrix) {
      if (!params.use_norm_as_weight)
        TRIQS_RUNTIME_ERROR << "To measure the density_matrix of atomic states, you need to set "
                               "use_norm_as_weight to True, i.e. to reweight the QMC";
      qmc.add_measure(measure_density_matrix{data, _density_matrix},
                      "Density Matrix for local static observable");
    }

    qmc.add_measure(measure_average_sign{data, _average_sign}, "Average sign");
    qmc.add_measure(measure_average_order{data, _average_order}, "Average order");
    qmc.add_measure(measure_auto_corr_time{data, _auto_corr_time}, "Auto-correlation time");

    // --------------------------------------------------------------------------

    // Run! The empty (starting) configuration has sign = 1
    _solve_status =
       qmc.warmup_and_accumulate(params.n_warmup_cycles, params.n_cycles, params.length_cycle,
                                 triqs::utility::clock_callback(params.max_time));
    qmc.collect_results(_comm);

    if (params.verbosity >= 2) {
      std::cout << "Average sign: " << _average_sign << std::endl;
      std::cout << "Average order: " << _average_order << std::endl;
      std::cout << "Auto-correlation time: " << _auto_corr_time << std::endl;
    }

    // Copy local (real or complex) G_tau back to complex G_tau
    if (G_tau && G_tau_accum) *G_tau = *G_tau_accum;
  }
} // namespace triqs_cthyb
