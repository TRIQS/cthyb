#pragma once
namespace cthyb {

 using namespace triqs::utility;
 using real_operator_t = many_body_operator<double>;

// All the arguments of the solve function
struct solve_parameters_t {

 /// Atomic Hamiltonian
 real_operator_t h_loc;

 /// Number of QMC cycles
 int n_cycles;

 /// Quantum numbers
 std::vector<real_operator_t> quantum_numbers = std::vector<real_operator_t>{};

 /// Length of a single QMC cycle
 int length_cycle = 50;

 /// Number of cycles for thermalization
 int n_warmup_cycles = 5000;

 /// Seed for random number generator
 int random_seed = 34788 + 928374 * boost::mpi::communicator().rank();

 /// Name of random number generator
 std::string random_name = "";

 /// Maximum runtime in seconds, use -1 to set infinite
 int max_time = -1;

 /// Verbosity level
 int verbosity = ((boost::mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes

 /// Calculate the full trace or use an estimate?
 bool use_trace_estimator = false;

 /// Whether to measure G(tau)
 bool measure_g_tau = true;

 /// Whether to measure G_l (Legendre)
 bool measure_g_l = false;

 /// Whether to measure perturbation order
 bool measure_pert_order = false;

 /// Make the analysis histograms of the trace computation
 bool make_histograms = false;

 solve_parameters_t() {}
 
 solve_parameters_t(real_operator_t h_loc, int n_cycles) : h_loc(h_loc), n_cycles(n_cycles) {}

};
}
