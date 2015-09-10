#pragma once

namespace cthyb {

 using namespace triqs::operators;
 using real_operator_t = many_body_operator;

// All the arguments of the solve function
struct solve_parameters_t {

 /// Interacting part of the atomic Hamiltonian
 /// type: Operator
 real_operator_t h_int;

 /// Number of QMC cycles
 int n_cycles;

 /// Partition method
 /// type: str
 std::string partition_method = "autopartition";

 /// Quantum numbers
 /// type: list(Operator)
 /// default: []
 std::vector<real_operator_t> quantum_numbers = std::vector<real_operator_t>{};

 /// Length of a single QMC cycle
 int length_cycle = 50;

 /// Number of cycles for thermalization
 int n_warmup_cycles = 5000;

 /// Seed for random number generator
 /// default: 34788 + 928374 * MPI.rank
 int random_seed = 34788 + 928374 * triqs::mpi::communicator().rank();

 /// Name of random number generator
 /// type: str
 std::string random_name = "";

 /// Maximum runtime in seconds, use -1 to set infinite
 int max_time = -1;

 /// Verbosity level
 /// default: 3 on MPI rank 0, 0 otherwise.
 int verbosity = ((triqs::mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes

 /// Add shifting a move as a move?
 bool move_shift = true;

 /// Add double insertions as a move?
 bool move_double = false;

 /// Calculate the full trace or use an estimate?
 bool use_trace_estimator = false;

 /// Measure G(tau)?
 bool measure_g_tau = true;

 /// Measure G_l (Legendre)?
 bool measure_g_l = false;

 /// Measure perturbation order?
 bool measure_pert_order = false;

 /// Measure the contribution of each atomic state to the trace?
 bool measure_density_matrix= false;

 /// Use the norm of the density matrix in the weight if true, otherwise use Trace
 bool use_norm_as_weight = false;

 /// Analyse performance of trace computation with histograms (developers only)?
 bool performance_analysis = false;

 /// Operator insertion/removal probabilities for different blocks
 /// type: dict(str:float)
 /// default: {}
 std::map<std::string,double> proposal_prob = (std::map<std::string,double>{});

 solve_parameters_t() {}

 solve_parameters_t(real_operator_t h_int, int n_cycles) : h_int(h_int), n_cycles(n_cycles) {}

};
}
