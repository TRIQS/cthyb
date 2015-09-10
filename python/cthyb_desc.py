# Generated automatically using the command :
# c++2py.py ../c++/solver_core.hpp -p -mpytriqs.applications.impurity_solvers.cthyb -o cthyb --moduledoc "The cthyb solver"
from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.applications.impurity_solvers.cthyb", doc = "The cthyb solver")

# All the triqs C++/Python modules
module.use_module('gf')
module.use_module('operators')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/solver_core.hpp")
module.add_include("../c++/post_process.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/pair.hpp>
#include <triqs/python_tools/converters/map.hpp>
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/tuple.hpp>
using namespace triqs::gfs;
using triqs::operators::many_body_operator;
using namespace cthyb;
#include "./converters.hxx"
""")

# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "solver_core",   # name of the C++ class
)

c.add_constructor("""(double beta, std::map<std::string,indices_type> gf_struct, int n_iw = 1025, int n_tau = 10001, int n_l = 50)""",
                  doc = """ """)

c.add_method("""void solve (**cthyb::solve_parameters_t)""",
             doc = """  Parameter Name                            Type            Default                       Documentation

  h_int                                     Operator        --                            Interacting part of the atomic Hamiltonian
  n_cycles                                  int             --                            Number of QMC cycles
  partition_method                          str             "autopartition"               Partition method
  quantum_numbers                           list(Operator)  []                            Quantum numbers
  length_cycle                              int             50                            Length of a single QMC cycle
  n_warmup_cycles                           int             5000                          Number of cycles for thermalization
  random_seed                               int             34788 + 928374 * MPI.rank     Seed for random number generator
  random_name                               str             ""                            Name of random number generator
  max_time                                  int             -1                            Maximum runtime in seconds, use -1 to set infinite
  verbosity                                 int             3 on MPI rank 0, 0 otherwise. Verbosity level
  move_shift                                bool            true                          Add shifting a move as a move?
  move_double                               bool            false                         Add double insertions as a move?
  use_trace_estimator                       bool            false                         Calculate the full trace or use an estimate?
  measure_g_tau                             bool            true                          Measure G(tau)?
  measure_g_l                               bool            false                         Measure G_l (Legendre)?
  measure_pert_order                        bool            false                         Measure perturbation order?
  measure_density_matrix                    bool            false                         Measure the contribution of each atomic state to the trace?
  use_norm_as_weight_for_atomic_correlators bool            false                         Use the norm of the density matrix in the weight if true, otherwise use Trace
  performance_analysis                      bool            false                         Analyse performance of trace computation with histograms (developers only)?
  proposal_prob                             dict(str:float) {}                            Operator insertion/removal probabilities for different blocks                  """)

#c.add_method("""std::map<std::string,std::vector<double>> atomic_observables (std::map<std::string,real_operator_t> obs_map)""",
#             doc = """Static observables of the atomic problem """)

c.add_property(name = "h_loc",
               getter = cfunction("many_body_op_type h_loc ()"),
               doc = """The local Hamiltonian of the problem """)

c.add_property(name = "last_solve_parameters",
               getter = cfunction("cthyb::solve_parameters_t get_last_solve_parameters ()"),
               doc = """Set of parameters used in the last call to solve """)

c.add_property(name = "G0_iw",
               getter = cfunction("block_gf_view<imfreq> G0_iw ()"),
               doc = """G0(iw) in imaginary frequencies """)

c.add_property(name = "Delta_tau",
               getter = cfunction("block_gf_view<imtime> Delta_tau ()"),
               doc = """Delta(tau) in imaginary time """)

c.add_property(name = "G_tau",
               getter = cfunction("block_gf_view<imtime> G_tau ()"),
               doc = """G(tau) in imaginary time """)

c.add_property(name = "G_l",
               getter = cfunction("block_gf_view<legendre> G_l ()"),
               doc = """G_l in Legendre polynomials representation """)

c.add_property(name = "atomic_gf",
               getter = cfunction("block_gf_view<imtime> atomic_gf ()"),
               doc = """Atomic G(tau) in imaginary time """)

c.add_property(name = "density_matrix",
               getter = cfunction("std::vector<matrix<double>> density_matrix ()"),
               doc = """Density matrix """)

c.add_method("""h_loc_diagonalization get_h_loc_diagonalization ()""",
             doc = """Access to the Hloc diagonalization """)

c.add_property(name = "average_sign",
               getter = cfunction("mc_sign_type average_sign ()"),
               doc = """Monte Carlo average sign """)

c.add_property(name = "solve_status",
               getter = cfunction("int solve_status ()"),
               doc = """Status of the solve on exit """)

module.add_class(c)

#The class sorted_spaces
c = class_(
        py_type = "SortedSpaces",  # name of the python class
        c_type = "sorted_spaces",   # name of the C++ class
        hdf5 = True,
)

c.add_constructor("""(many_body_op_t h_, std::vector<many_body_op_t> qn_vector)""",
                  doc = """ """)

c.add_constructor("""(many_body_op_t h_, std::vector<indices_t> gf_struct)""",
                  doc = """ """)

c.add_constructor("""(many_body_op_t h_, std::vector<many_body_op_t> qn_vector, std::vector<indices_t> gf_struct)""",
                  doc = """ """)

c.add_method("""block_matrix_t matrix_element (std::vector<int> P)""",
             doc = """For a symmetry implemented as a permutation of the C operators, return the matrix by block
 Blocks must be invariant.
 The permutation is a vector i-> j, where the int are the indices of the fundamental_operator_set of the object. """)

c.add_method("""block_matrix_t matrix_element (std::vector<std::pair<indices_t,indices_t>> P)""",
             doc = """ """)

c.add_method("""int dim ()""",
             doc = """Dimension of the full Hilbert space """)

c.add_method("""int n_blocks ()""",
             doc = """Number of Blocks """)

c.add_method("""int get_block_dim (int b)""",
             doc = """The dimension of block b """)

c.add_method("""int n_c_operators ()""",
             doc = """Number of c operators """)

c.add_method("""long fundamental_operator_connect (bool dagger, int index, int n)""",
             doc = """Connections for fundamental operators """)

c.add_method("""matrix<double> fundamental_operator_matrix (bool dagger, int index, int block_index)""",
             doc = """Connections for fundamental operators """)

c.add_method("""double get_eigenvalue (int bl, int i)""",
             doc = """Get the i-th eigenvalue of block bl """)

c.add_method("""int flatten_block_index (int bl, int i)""",
             doc = """bl : block index, i : index within the block. Returns the index in the full hilbert space """)

c.add_method("""int get_vacuum_block_index ()""",
             doc = """Vaccum is necessarly a block of size 1. Returns the block index. """)

c.add_method("""full_hilbert_space_state_t get_vacuum ()""",
             doc = """Vaccum is necessarly a block of size 1. Returns the block index. """)

c.add_method("""double average (block_matrix_t density_matrix, many_body_op_t op)""",
             doc = """Trace (op * density_matrix) """)

c.add_method("""full_hilbert_space_state_t act (many_body_op_t op, full_hilbert_space_state_t st)""",
             doc = """Act with operator Op on state St """)

c.add_method("""double average_on_projector (block_matrix_t density_matrix, full_hilbert_space_state_t psi)""",
             doc = """Average the density matrix on the state psi """)

c.add_method("""std::string eigenstate_repr (int bl, int k)""",
             doc = """A human readable representation for eigenstate """)

c.add_method("""double get_gs_energy ()""",
             doc = """Ground state energy (i.e. min of all subspaces). """)

c.add_method("""std::vector<std::tuple<vector<double>,std::vector<double>>> get_energy_and_quantum_numbers ()""",
             doc = """ """)

c.add_method("""double partition_function (double beta)""",
             doc = """The partition function """)

c.add_method("""block_gf<imtime> atomic_gf (double beta, std::map<std::string,many_body_op_t::indices_t> gf_struct, int n_tau, std::vector<std::pair<int, int>>  excluded_states)""",
             doc = """The atomic green function """)

module.add_class(c)

module.generate_code()
