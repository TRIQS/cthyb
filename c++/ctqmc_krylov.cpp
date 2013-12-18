#include "ctqmc_krylov.hpp"
#include <triqs/utility/callbacks.hpp>

#include "move_insert.hpp"
#include "move_remove.hpp"
#include "measure_g.hpp"

namespace cthyb_krylov {

void ctqmc_krylov::solve(utility::parameters p_in) {

 p_in.update(solve_defaults());
 auto const& params = p_in;

 qmc_data data(params, sosp, deltat);
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
 
 // run!! The empty configuration has sign = 1
 qmc.start(1.0, triqs::utility::clock_callback(params["max_time"]));
 qmc.collect_results(c);
}

//----------------------------------------------------------------------------------------------
parameter_defaults ctqmc_krylov::constructor_defaults() const {

 parameter_defaults pdef;

 pdef.required("beta", double(), "Inverse temperature")
     .optional("n_tau_delta", int(10001), "Number of time slices for Delta(tau)")
     .optional("n_tau_g", int(10001), "Number of time slices for G(tau)")
     .optional("n_w", int(1025), "Number of Matsubara frequencies");

 return pdef;
}

//----------------------------------------------------------------------------------------------

parameter_defaults ctqmc_krylov::solve_defaults() const {

 parameter_defaults pdef;

 pdef.required("n_cycles", int(), "Number of QMC cycles")
     .optional("length_cycle", int(50), "Length of a single QMC cycle")
     .optional("n_warmup_cycles", int(5000), "Number of cycles for thermalization")
     .optional("random_seed", int(34788 + 928374 * c.rank()), "Seed for random number generator")
     .optional("random_name", std::string(""), "Name of random number generator")
     .optional("max_time", int(-1), "Maximum runtime in seconds, use -1 to set infinite")
     .optional("verbosity", int(3), "Verbosity level")
     .optional("measure_gt", bool(true), "Whether to measure G(tau)")
     .optional("make_path_histograms", bool(false), " Make the analysis histograms of the trace computation ")
     .optional("krylov_bs_use_cutoff", bool(false), " bool ")
     .optional("krylov_bs_prob_cutoff", double(1e-8), " double ") // put negative to include all boundary states.
     .optional("krylov_gs_energy_convergence", 1e-10, " double ")
     .optional("krylov_small_matrix_size", int(10), " unsigned int ");
 return pdef;
}

void ctqmc_krylov::help() const {
 // TODO
}
}
