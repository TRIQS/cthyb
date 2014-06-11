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
#include "ctqmc.hpp"
#include <triqs/utility/callbacks.hpp>

#include "move_insert.hpp"
#include "move_remove.hpp"
#include "measure_g.hpp"
#include "measure_perturbation_hist.hpp"



namespace cthyb {

ctqmc::ctqmc(double beta_, std::map<std::string,std::vector<int>> const & gf_struct_, int n_tau_delta, int n_tau_g):
  beta(beta_), gf_struct(gf_struct_) {

  std::vector<std::string> block_names;
  std::vector<gf<imtime>> deltat_blocks;
  std::vector<gf<imtime>> gt_blocks;

  for (auto const& block : gf_struct) {
    block_names.push_back(block.first);
    int n = block.second.size();
    deltat_blocks.push_back(gf<imtime>{{beta, Fermion, n_tau_delta, half_bins}, {n, n}});
    gt_blocks.push_back(gf<imtime>{{beta, Fermion, n_tau_g, half_bins}, {n, n}});
  }

  deltat = make_block_gf(block_names, deltat_blocks);
  gt = make_block_gf(block_names, gt_blocks);

}


void ctqmc::solve(real_operator_t const & h_loc, params::parameters params,
                  std::vector<real_operator_t> const & quantum_numbers, bool use_quantum_numbers) {

  // basis of operators to use
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

  if (use_quantum_numbers)
   sosp = {h_loc, quantum_numbers, fops};
  else 
   sosp = {h_loc, fops};

 qmc_data data(beta, params, sosp, linindex, deltat);
 mc_tools::mc_generic<mc_sign_type> qmc(params);

 // Moves
 auto& delta_names = deltat.domain().names();
 for (size_t block = 0; block < deltat.domain().size(); ++block) {
  int block_size = deltat[block].data().shape()[1];
  qmc.add_move(move_insert_c_cdag(block, block_size, data, qmc.rng(), false), "Insert Delta_" + delta_names[block]);
  qmc.add_move(move_remove_c_cdag(block, block_size, data, qmc.rng()), "Remove Delta_" + delta_names[block]);
 }

 // Measurements
 if (params["measure_gt"]) {
  auto& gt_names = gt.domain().names();
  for (size_t block = 0; block < gt.domain().size(); ++block) {
   qmc.add_measure(measure_g(block, gt[block], data), "G measure (" + gt_names[block] + ")");
  }
 }
 if (params["measure_pert_order"]) {
  auto& gt_names = gt.domain().names();
  for (size_t block = 0; block < gt.domain().size(); ++block) {
   qmc.add_measure(measure_perturbation_hist(block, data, "histo_pert_order_" + gt_names[block] + ".dat"), "Perturbation order (" + gt_names[block] + ")");
  }
 }

 // run!! The empty configuration has sign = 1
 qmc.start(1.0, triqs::utility::clock_callback(params["max_time"]));
 qmc.collect_results(c);
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
     .add_field("measure_gt", bool(true), "Whether to measure G(tau)")
     .add_field("measure_pert_order", bool(false), "Whether to measure perturbation order")
     .add_field("make_histograms", bool(false), "Make the analysis histograms of the trace computation");

 return pdef;
}

//----------------------------------------------------------------------------------------------

void ctqmc::help() {
 // TODO
}

}
