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
#include <triqs/utility/callbacks.hpp>
#include <triqs/utility/exceptions.hpp>
#include <triqs/utility/variant_int_string.hpp>
#include <triqs/gfs.hpp>
#include <fstream>

#include "move_insert.hpp"
#include "move_remove.hpp"
#include "move_double_insert.hpp"
#include "move_double_remove.hpp"
#include "move_shift.hpp"
#include "measure_g.hpp"
#include "measure_g_legendre.hpp"
#include "measure_perturbation_hist.hpp"
#include "measure_density_matrix.hpp"
#include "measure_average_sign.hpp"

namespace cthyb {

struct index_visitor  {
  std::vector<std::string> indices;
  void operator()(int i) { indices.push_back(std::to_string(i)); }
  void operator()(std::string s) { indices.push_back(s); }
};

solver_core::solver_core(double beta_, std::map<std::string, indices_type> const & gf_struct_, int n_iw, int n_tau, int n_l):
  beta(beta_), gf_struct(gf_struct_) {

  if ( n_tau < 2*n_iw ) TRIQS_RUNTIME_ERROR << "Must use as least twice as many tau points as Matsubara frequencies: n_iw = " << n_iw << " but n_tau = " << n_tau << ".";

  std::vector<std::string> block_names;
  std::vector<gf<imfreq>> g0_iw_blocks;
  std::vector<gf<imtime>> g_tau_blocks;
  std::vector<gf<legendre>> g_l_blocks;
  std::vector<gf<imtime>> delta_tau_blocks;
  std::vector<gf<imtime,matrix_real_valued>> g_tau_real_blocks; // Local real quantities for accumulation

  for (auto const& bl : gf_struct) {
    block_names.push_back(bl.first);
    int n = bl.second.size();

    index_visitor iv;
    for (auto & ind: bl.second) { apply_visitor(iv, ind); }
    std::vector<std::vector<std::string>> indices{{iv.indices,iv.indices}};

    g0_iw_blocks.push_back(gf<imfreq>{{beta, Fermion, n_iw}, {n, n}, indices});
    g_tau_blocks.push_back(gf<imtime>{{beta, Fermion, n_tau}, {n, n}, indices});
    g_l_blocks.push_back(gf<legendre>{{beta, Fermion, static_cast<size_t>(n_l)}, {n,n}, indices});
    delta_tau_blocks.push_back(gf<imtime>{{beta, Fermion, n_tau}, {n, n}, indices});
    g_tau_real_blocks.push_back(gf<imtime,matrix_real_valued>{{beta, Fermion, n_tau}, {n, n}});
  }

  _G0_iw = make_block_gf(block_names, g0_iw_blocks);
  _G_tau = make_block_gf(block_names, g_tau_blocks);
  _G_l = make_block_gf(block_names, g_l_blocks);
  _Delta_tau = make_block_gf(block_names, delta_tau_blocks);
  _G_tau_real = make_block_gf(block_names, g_tau_real_blocks);

}

/// -------------------------------------------------------------------------------------------

void solver_core::solve(solve_parameters_t const & params) {

  _last_solve_parameters = params;

  // determine basis of operators to use
  fundamental_operator_set fops;
  for (auto const & bl: gf_struct) {
    for (auto const & a: bl.second) {
      fops.insert(bl.first, a);
    }
  }

  // setup the linear index map
  std::map<std::pair<int,int>,int> linindex;
  int block_index = 0;
  for (auto const & bl: gf_struct) {
    int inner_index = 0;
    for (auto const & a: bl.second) {
      linindex[std::make_pair(block_index, inner_index)] = fops[{bl.first,a}];
      inner_index++;
    }
    block_index++;
  }

  // Make list of block sizes
  std::vector<int> n_inner;
  for (auto const& bl : gf_struct) {
    n_inner.push_back(bl.second.size());
  }

  // Calculate imfreq quantities
  auto G0_iw_inv = map([](gf_const_view<imfreq> x){return triqs::gfs::inverse(x);}, _G0_iw);
  auto Delta_iw = G0_iw_inv;

  _h_loc = params.h_int;

  // Add quadratic terms to h_loc
  int b = 0;
  for (auto const & bl: gf_struct) {
    int n1 = 0;
    for (auto const & a1: bl.second) {
      int n2 = 0;
      for (auto const & a2: bl.second) {
        _h_loc = _h_loc + _G0_iw[b].singularity()(2)(n1,n2).real() * c_dag<double>(bl.first,a1)*c<double>(bl.first,a2);
        n2++;
      }
      n1++;
    }
    b++;
  }

  // Determine terms Delta_iw from G0_iw and ensure that the 1/iw behaviour of G0_iw is correct
  b = 0;
  triqs::clef::placeholder<0> iw_;
  for (auto const & bl: gf_struct) {
    Delta_iw[b](iw_) << G0_iw_inv[b].singularity()(-1)*iw_ + G0_iw_inv[b].singularity()(0);
    Delta_iw[b] = Delta_iw[b] - G0_iw_inv[b];
    _Delta_tau[b]() = inverse_fourier(Delta_iw[b]);
    b++;
  }

  // Report what h_loc we are using
  if (params.verbosity >= 2) std::cout << "The local Hamiltonian of the problem:" << std::endl << _h_loc << std::endl;

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
   if(params.measure_density_matrix) _density_matrix = atomic_density_matrix(h_diag, beta);
   return;
  }

  // Initialise Monte Carlo quantities
  qmc_data data(beta, params, h_diag, linindex, _Delta_tau, n_inner);
  auto qmc = mc_tools::mc_generic<mc_sign_type>(params.n_cycles, params.length_cycle, params.n_warmup_cycles, params.random_name,
                                                params.random_seed, params.verbosity);

  // Moves
  using move_set_type = mc_tools::move_set<mc_sign_type>;
  move_set_type inserts(qmc.get_rng());
  move_set_type removes(qmc.get_rng());
  move_set_type double_inserts(qmc.get_rng());
  move_set_type double_removes(qmc.get_rng());

  auto& delta_names = _Delta_tau.domain().names();
  for (size_t block = 0; block < _Delta_tau.domain().size(); ++block) {
   int block_size = _Delta_tau[block].data().shape()[1];
   auto const& block_name = delta_names[block];
   auto f = params.proposal_prob.find(block_name);
   double prop_prob = (f != params.proposal_prob.end() ? f->second : 1.0);
   inserts.add(move_insert_c_cdag(block, block_size, data, qmc.get_rng(), params.performance_analysis), "Insert Delta_" + block_name, prop_prob);
   removes.add(move_remove_c_cdag(block, block_size, data, qmc.get_rng(), params.performance_analysis), "Remove Delta_" + block_name, prop_prob);
   if (params.move_double) {
    for (size_t block2 = 0; block2 < _Delta_tau.domain().size(); ++block2) {
     int block_size2 = _Delta_tau[block2].data().shape()[1];
     double_inserts.add(move_insert_c_c_cdag_cdag(block, block2, block_size, block_size2, data, qmc.get_rng(), params.performance_analysis),
                 "Insert Delta_" + delta_names[block] + "_" + delta_names[block2], 1.0);
     double_removes.add(move_remove_c_c_cdag_cdag(block, block2, block_size, block_size2, data, qmc.get_rng(), params.performance_analysis),
                 "Remove Delta_" + delta_names[block] + "_" + delta_names[block2], 1.0);
    }
   }
  }

  qmc.add_move(inserts, "Insert two operators", 1.0);
  qmc.add_move(removes, "Remove two operators", 1.0);
  if (params.move_double) {
   qmc.add_move(double_inserts, "Insert four operators", 1.0);
   qmc.add_move(double_removes, "Remove four operators", 1.0);
  }
  if (params.move_shift) qmc.add_move(move_shift_operator(data, qmc.get_rng(), params.performance_analysis), "Shift one operator", 1.0);
 
  // Measurements
  if (params.measure_g_tau) {
   auto& g_names = _G_tau.domain().names();
   for (size_t block = 0; block < _G_tau.domain().size(); ++block) {
    qmc.add_measure(measure_g(block, _G_tau_real[block], data), "G measure (" + g_names[block] + ")");
   }
  }
  if (params.measure_g_l) {
   auto& g_names = _G_l.domain().names();
   for (size_t block = 0; block < _G_l.domain().size(); ++block) {
    qmc.add_measure(measure_g_legendre(block, _G_l[block], data), "G_l measure (" + g_names[block] + ")");
   }
  }
  if (params.measure_pert_order) {
   auto& g_names = _G_tau.domain().names();
   for (size_t block = 0; block < _G_tau.domain().size(); ++block) {
    qmc.add_measure(measure_perturbation_hist(block, data, "histo_pert_order_" + g_names[block] + ".dat"), "Perturbation order (" + g_names[block] + ")");
   }
   qmc.add_measure(measure_perturbation_hist_total(data, "histo_pert_order.dat"), "Perturbation order");
  }

  if (params.measure_density_matrix) {
   if (!params.use_norm_as_weight)
    TRIQS_RUNTIME_ERROR << "To measure the density_matrix of atomic states, you need to set "
                           "use_norm_as_weight to True, i.e. to reweight the QMC";
   qmc.add_measure(measure_density_matrix{data, _density_matrix}, "Density Matrix for local static observable");
  }

  qmc.add_measure(measure_average_sign{data, _average_sign}, "Average sign");

  // Run! The empty (starting) configuration has sign = 1
  _solve_status = qmc.start(1.0, triqs::utility::clock_callback(params.max_time));
  qmc.collect_results(_comm);

  std::cout << "Average sign: " << _average_sign << std::endl;

  // Copy real G_tau back to complex G_tau
  if (params.measure_g_tau) {
   _G_tau = _G_tau_real;
  }

}

}
