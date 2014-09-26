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
#include "ctqmc.hpp"
#include <triqs/utility/callbacks.hpp>
#include <triqs/utility/exceptions.hpp>
#include <triqs/gfs.hpp>

#include "move_insert.hpp"
#include "move_remove.hpp"
#include "measure_g.hpp"
#include "measure_g_legendre.hpp"
#include "measure_perturbation_hist.hpp"

namespace cthyb {

ctqmc::ctqmc(double beta_, std::map<std::string,std::vector<int>> const & gf_struct_, int n_iw, int n_tau, int n_l):
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

void ctqmc::solve(real_operator_t h_loc, params::parameters params,
                  std::vector<real_operator_t> const & quantum_numbers, bool use_quantum_numbers) {

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

  // Calculate imfreq quantities
  auto G0_iw_inv = map([](gf_const_view<imfreq> x){return triqs::gfs::inverse(x);}, _G0_iw);
  auto Delta_iw = G0_iw_inv;

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

  // Determine block structure
  if (use_quantum_numbers)
   sosp = {h_loc, quantum_numbers, fops};
  else 
   sosp = {h_loc, fops};

  if(params["make_histograms"]) std::ofstream("Diagonalization_atomic_pb") << sosp;

  qmc_data data(beta, params, sosp, linindex, _Delta_tau);
  mc_tools::mc_generic<mc_sign_type> qmc(params);
 
  // Moves
  auto& delta_names = _Delta_tau.domain().names();
  for (size_t block = 0; block < _Delta_tau.domain().size(); ++block) {
   int block_size = _Delta_tau[block].data().shape()[1];
   qmc.add_move(move_insert_c_cdag(block, block_size, data, qmc.rng(), false), "Insert Delta_" + delta_names[block]);
   qmc.add_move(move_remove_c_cdag(block, block_size, data, qmc.rng()), "Remove Delta_" + delta_names[block]);
  }
 
  // Measurements
  if (params["measure_g_tau"]) {
   auto& g_names = _G_tau.domain().names();
   for (size_t block = 0; block < _G_tau.domain().size(); ++block) {
    qmc.add_measure(measure_g(block, _G_tau[block], data), "G measure (" + g_names[block] + ")");
   }
  }
  if (params["measure_g_l"]) {
   auto& g_names = _G_l.domain().names();
   for (size_t block = 0; block < _G_l.domain().size(); ++block) {
    qmc.add_measure(measure_g_legendre(block, _G_l[block], data), "G_l measure (" + g_names[block] + ")");
   }
  }
  if (params["measure_pert_order"]) {
   auto& g_names = _G_tau.domain().names();
   for (size_t block = 0; block < _G_tau.domain().size(); ++block) {
    qmc.add_measure(measure_perturbation_hist(block, data, "histo_pert_order_" + g_names[block] + ".dat"), "Perturbation order (" + g_names[block] + ")");
   }
  }
 
  // run!! The empty configuration has sign = 1
  qmc.start(1.0, triqs::utility::clock_callback(params["max_time"]));
  qmc.collect_results(_comm);
}

//----------------------------------------------------------------------------------------------

parameters_t ctqmc::solve_parameters() {

 auto pdef = parameters_t{};
 boost::mpi::communicator world;

 pdef.add_field("n_cycles", no_default<int>(), "Number of QMC cycles")
     .add_field("length_cycle", int(50), "Length of a single QMC cycle")
     .add_field("n_warmup_cycles", int(5000), "Number of cycles for thermalization")
     .add_field("random_seed", int(34788 + 928374 * world.rank()), "Seed for random number generator")
     .add_field("random_name", std::string(""), "Name of random number generator")
     .add_field("max_time", int(-1), "Maximum runtime in seconds, use -1 to set infinite")
     .add_field("verbosity", (world.rank()==0 ? int(3) : int(0)), "Verbosity level")
     .add_field("use_trace_estimator", bool(false), "Calculate the full trace or use an estimate?")
     .add_field("measure_g_tau", bool(true), "Whether to measure G(tau)")
     .add_field("measure_g_l", bool(false), "Whether to measure G_l (Legendre)")
     .add_field("measure_pert_order", bool(false), "Whether to measure perturbation order")
     .add_field("make_histograms", bool(false), "Make the analysis histograms of the trace computation");

 return pdef;
}

//----------------------------------------------------------------------------------------------

void ctqmc::help() {
 // TODO
}

using namespace triqs::gfs;

gf<imtime> change_mesh(gf_const_view<imtime> old_gf, int new_n_tau) {
    auto const& old_m = old_gf.mesh();
    gf<imtime> new_gf{{old_m.domain().beta, old_m.domain().statistic, new_n_tau, old_m.kind()}, old_gf.get_data_shape().front_pop()};
    auto const& new_m = new_gf.mesh();

    new_gf.data()() = 0;
    double f = old_m.delta()/new_m.delta();
    for(auto tau : old_m) new_gf[closest_mesh_pt(double(tau))] += f*old_gf[tau];

    new_gf.singularity() = old_gf.singularity();

    return new_gf;
}

block_gf<imtime> change_mesh(block_gf_const_view<imtime> old_gf, int new_n_tau) {
    std::vector<gf<imtime>> blocks;
    for(auto const& bl : old_gf) blocks.push_back(change_mesh(bl,new_n_tau));
    return make_block_gf(old_gf.domain().names(), blocks);
}

}
