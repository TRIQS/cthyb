/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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
#include <fstream>
#include "solver_core.hpp"
#include <triqs/utility/callbacks.hpp>
#include <triqs/utility/exceptions.hpp>
#include <triqs/gfs.hpp>

#include "move_insert.hpp"
#include "move_remove.hpp"
#include "move_shift.hpp"
#include "measure_g.hpp"
#include "measure_g_legendre.hpp"
#include "measure_perturbation_hist.hpp"

namespace cthyb {

solver_core::solver_core(double beta_, std::map<std::string,std::vector<int>> const & gf_struct_, int n_iw, int n_tau, int n_l):
  beta(beta_), gf_struct(gf_struct_) {

  if ( n_tau < 2*n_iw ) TRIQS_RUNTIME_ERROR << "Must use as least twice as many tau points as Matsubara frequencies: n_iw = " << n_iw << " but n_tau = " << n_tau << ".";

  std::vector<std::string> block_names;
  std::vector<gf<imfreq>> g0_iw_blocks;
  std::vector<gf<imtime>> g_tau_blocks;
  std::vector<gf<legendre>> g_l_blocks;
  std::vector<gf<imtime>> delta_tau_blocks;

  for (auto const& block : gf_struct) {
    block_names.push_back(block.first);
    int n = block.second.size();
    g0_iw_blocks.push_back(gf<imfreq>{{beta, Fermion, n_iw}, {n, n}});
    g_tau_blocks.push_back(gf<imtime>{{beta, Fermion, n_tau, full_bins}, {n, n}});
    g_l_blocks.push_back(gf<legendre>{{beta, Fermion, static_cast<size_t>(n_l)}, {n,n}});
    delta_tau_blocks.push_back(gf<imtime>{{beta, Fermion, n_tau, full_bins}, {n, n}});
  }

  _G0_iw = make_block_gf(block_names, g0_iw_blocks);
  _G_tau = make_block_gf(block_names, g_tau_blocks);
  _G_l = make_block_gf(block_names, g_l_blocks);
  _Delta_tau = make_block_gf(block_names, delta_tau_blocks);

}

// TODO Move functions below to triqs library
template<typename F, typename T>
std::vector<std14::result_of_t<F(T)>> map(F && f, std::vector<T> const & V) {
  std::vector<std14::result_of_t<F(T)>> res;
  res.reserve(V.size());
  for(auto & x : V) res.emplace_back(f(x));
  return res;
}

template<typename F, typename G>
gf<block_index,std14::result_of_t<F(G)>> map(F && f, gf<block_index,G> const & g) {
  return make_block_gf(map(f, g.data()));
}

/// -------------------------------------------------------------------------------------------

void solver_core::solve(solve_parameters_t const & params) { 

  // determine basis of operators to use
  fundamental_operator_set fops;
  std::map<std::pair<int,int>,int> linindex;
  int block_index = 0;
  for (auto const & b: gf_struct) {
    int inner_index = 0;
    for (auto const & a: b.second) {
      fops.insert(b.first, a);
      linindex[std::make_pair(block_index, inner_index)] = fops[{b.first,a}];
      inner_index++;
    }
    block_index++;
  }

  // Make list of block sizes
  std::vector<int> n_inner;
  for (auto const& block : gf_struct) {
    n_inner.push_back(block.second.size());
  }

  // Calculate imfreq quantities
  auto G0_iw_inv = map([](gf_const_view<imfreq> x){return triqs::gfs::inverse(x);}, _G0_iw);
  auto Delta_iw = G0_iw_inv;
  auto h_loc = params.h_loc;

  // Add quadratic terms to h_loc
  int b = 0;
  for (auto const & bl: gf_struct) {
    for (auto const & a1: bl.second) {
      for (auto const & a2: bl.second) {
        h_loc = h_loc + _G0_iw[b].singularity()(2)(a1,a2).real() * c_dag(bl.first,a1)*c(bl.first,a2);
      }
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
    _G0_iw[b](iw_) << iw_ + G0_iw_inv[b].singularity()(0) ;
    _G0_iw[b] = _G0_iw[b] - Delta_iw[b];
    _G0_iw[b]() = triqs::gfs::inverse(_G0_iw[b]);
    b++;
  }

  // Report what h_loc we are using
  if (params.verbosity >= 2) std::cout << "The local Hamiltonian of the problem:" << std::endl << h_loc << std::endl;

  // Determine block structure
  if (params.partition_method == "autopartition") {
   if (params.verbosity >= 2) std::cout << "Using autopartition algorithm to partition the local Hilbert space" << std::endl;
   sosp = {h_loc, fops};
  } else if (params.partition_method == "quantum_numbers") {
   if (params.quantum_numbers.empty()) TRIQS_RUNTIME_ERROR << "No quantum numbers provided.";
   if (params.verbosity >= 2) std::cout << "Using quantum numbers to partition the local Hilbert space" << std::endl;
   sosp = {h_loc, params.quantum_numbers, fops};
  } else if (params.partition_method == "none") { // give empty quantum numbers list
   std::cout << "Not partitioning the local Hilbert space" << std::endl;
   sosp = {h_loc, std::vector<real_operator_t>{}, fops};
  } else
   TRIQS_RUNTIME_ERROR << "Partition method " << params.partition_method << " not recognised.";

  if (params.make_histograms) std::ofstream("Diagonalization_atomic_pb") << sosp;

  qmc_data data(beta, params, sosp, linindex, _Delta_tau, n_inner);
  auto qmc = mc_tools::mc_generic<mc_sign_type>(params.n_cycles, params.length_cycle, params.n_warmup_cycles, params.random_name,
                                                params.random_seed, params.verbosity);

  // Moves
  auto& delta_names = _Delta_tau.domain().names();
  for (size_t block = 0; block < _Delta_tau.domain().size(); ++block) {
   int block_size = _Delta_tau[block].data().shape()[1];
   qmc.add_move(move_insert_c_cdag(block, block_size, data, qmc.rng(), params.make_histograms), "Insert Delta_" + delta_names[block]);
   qmc.add_move(move_remove_c_cdag(block, block_size, data, qmc.rng()), "Remove Delta_" + delta_names[block]);
  }
  qmc.add_move(move_shift_operator(data, qmc.rng(), params.make_histograms), "Shift Operator");
  std::cout << " finished adding moves " << std::endl; // FIXME DEBUG
 
  // Measurements
  if (params.measure_g_tau) {
   auto& g_names = _G_tau.domain().names();
   for (size_t block = 0; block < _G_tau.domain().size(); ++block) {
    qmc.add_measure(measure_g(block, _G_tau[block], data), "G measure (" + g_names[block] + ")");
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
  }

  // Run! The empty configuration has sign = 1
  qmc.start(1.0, triqs::utility::clock_callback(params.max_time));
  qmc.collect_results(_comm);

}

//----------------------------------------------------------------------------------------------

using namespace triqs::gfs;

gf<imtime> change_mesh(gf_const_view<imtime> old_gf, int new_n_tau) {
    auto const& old_m = old_gf.mesh();
    gf<imtime> new_gf{{old_m.domain().beta, old_m.domain().statistic, new_n_tau, old_m.kind()}, get_target_shape(old_gf)};
    auto const& new_m = new_gf.mesh();

    new_gf.data()() = 0;
    double f = old_m.delta()/new_m.delta();
    for(auto tau : old_m) new_gf[closest_mesh_pt(double(tau))] += f*old_gf[tau];
    new_gf[0] *= 2.0; new_gf[new_n_tau-1] *= 2.0;

    new_gf.singularity() = old_gf.singularity();

    return new_gf;
}

block_gf<imtime> change_mesh(block_gf_const_view<imtime> old_gf, int new_n_tau) {
    std::vector<gf<imtime>> blocks;
    for(auto const& bl : old_gf) blocks.push_back(change_mesh(bl,new_n_tau));
    return make_block_gf(old_gf.domain().names(), blocks);
}

}
