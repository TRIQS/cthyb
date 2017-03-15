# Generated automatically using the command :
# c++2py.py ../c++/solver_core.hpp -I../../cthyb.build/c++ -I../c++ -p -mpytriqs.applications.impurity_solvers.cthyb -o cthyb --moduledoc "The cthyb solver"
from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.applications.impurity_solvers.cthyb", doc = "The cthyb solver")

# All the triqs C++/Python modules
module.use_module('gf', 'triqs')
module.use_module('multivar', 'triqs')
module.use_module('operators', 'triqs')
module.use_module('histograms', 'triqs')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("solver_core.hpp")
module.add_include("atom_diag.hpp")

module.add_enum("block_order", ["AABB","ABBA"], "cthyb", "Order of block indices for Block2Gf objects")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/gf.hpp>
#include <triqs/python_tools/converters/block_gf.hpp>
#include <triqs/python_tools/converters/block2_gf.hpp>
#include <triqs/python_tools/converters/pair.hpp>
#include <triqs/python_tools/converters/map.hpp>
#include <triqs/python_tools/converters/set.hpp>
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/variant.hpp>
#include <triqs/python_tools/converters/tuple.hpp>
#include <triqs/python_tools/converters/operators_real_complex.hpp>
#include <triqs/python_tools/converters/fundamental_operator_set.hpp>
using namespace triqs::gfs;
using triqs::operators::many_body_operator;
using namespace cthyb;
#include "./cthyb_converters.hxx"
""")

# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "solver_core",   # name of the C++ class
        doc = r"Core class of the cthyb solver",   # doc of the C++ class
)

c.add_constructor("""(double beta, std::map<std::string,indices_type> gf_struct, int n_iw = 1025, int n_tau = 10001, int n_l = 50)""",
                  doc = """ """)

solve_doc = """+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| Parameter Name         | Type                                           | Default                       | Documentation                                                                                           |
+========================+================================================+===============================+=========================================================================================================+
| h_int                  | Operator                                       |                               | Interacting part of the atomic Hamiltonian                                                              |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| n_cycles               | int                                            |                               | Number of QMC cycles                                                                                    |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| partition_method       | str                                            | "autopartition"               | Partition method                                                                                        |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| quantum_numbers        | list(Operator)                                 | []                            | Quantum numbers                                                                                         |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| length_cycle           | int                                            | 50                            | Length of a single QMC cycle                                                                            |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| n_warmup_cycles        | int                                            | 5000                          | Number of cycles for thermalization                                                                     |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| random_seed            | int                                            | 34788 + 928374 * MPI.rank     | Seed for random number generator                                                                        |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| random_name            | str                                            | ""                            | Name of random number generator                                                                         |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| max_time               | int                                            | -1 = infinite                 | Maximum runtime in seconds, use -1 to set infinite                                                      |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| verbosity              | int                                            | 3 on MPI rank 0, 0 otherwise. | Verbosity level                                                                                         |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| move_shift             | bool                                           | true                          | Add shifting an operator as a move?                                                                     |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| move_double            | bool                                           | false                         | Add double insertions as a move?                                                                        |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| use_trace_estimator    | bool                                           | false                         | Calculate the full trace or use an estimate?                                                            |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_g_tau          | bool                                           | true                          | Measure G(tau)?                                                                                         |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_g_l            | bool                                           | false                         | Measure G_l (Legendre)?                                                                                 |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_g2_inu         | bool                                           | false                         | Measure G^2(iomega,inu,inu')                                                                            |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_g2_legendre    | bool                                           | false                         | Measure G^2(iomega,l,l')                                                                                |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_g2_pp          | bool                                           | true                          | Measure G^2 within particle-particle channel.                                                           |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_g2_ph          | bool                                           | true                          | Measure G^2 within particle-hole channel.                                                               |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_g2_block_order | cthyb::block_order                             | AABB                          | Order of block indices in the definition of G^2.                                                        |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_g2_blocks      | std::set<std::pair<std::string, std::string> > | measure all blocks            | List of block index pairs of G^2 to measure.                                                            |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_g2_n_iw        | int                                            | 30                            | Number of bosonic Matsubara frequencies for G^2 measurement.                                            |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_g2_n_inu       | int                                            | 30                            | Number of fermionic Matsubara frequencies for G^2 measurement.                                          |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_g2_n_l         | int                                            | 20                            | Number of Legendre coefficients for G^2(iomega,l,l') measurement.                                       |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_pert_order     | bool                                           | false                         | Measure perturbation order?                                                                             |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| measure_density_matrix | bool                                           | false                         | Measure the reduced impurity density matrix?                                                            |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| use_norm_as_weight     | bool                                           | false                         | Use the norm of the density matrix in the weight if true, otherwise use Trace                           |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| performance_analysis   | bool                                           | false                         | Analyse performance of trace computation with histograms (developers only)?                             |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| proposal_prob          | dict(str:float)                                | {}                            | Operator insertion/removal probabilities for different blocks                                           |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| move_global            | dict(str : dict(indices : indices))            | {}                            | List of global moves (with their names). Each move is specified with an index substitution dictionary.  |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| move_global_prob       | double                                         | 0.05                          | Overall probability of the global moves                                                                 |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+
| imag_threshold         | double                                         | 1.e-15                        | Threshold below which imaginary components of Delta and h_loc are set to zero                           |
+------------------------+------------------------------------------------+-------------------------------+---------------------------------------------------------------------------------------------------------+ """

c.add_method("""void solve (**cthyb::solve_parameters_t)""", doc = solve_doc)

c.add_property(name = "h_loc",
               getter = cfunction("many_body_op_t h_loc ()"),
               doc = """The local Hamiltonian of the problem: :math:`H_{loc}` used in the last call to ``solve()``. """)

c.add_property(name = "last_solve_parameters",
               getter = cfunction("cthyb::solve_parameters_t last_solve_parameters ()"),
               doc = """Set of parameters used in the last call to ``solve()``. """)

c.add_property(name = "Delta_tau",
               getter = cfunction("block_gf_view<imtime> Delta_tau ()"),
               doc = """:math:`\\Delta(\\tau)` in imaginary time. """)

c.add_property(name = "G0_iw",
               getter = cfunction("block_gf_view<imfreq> G0_iw ()"),
               doc = """:math:`G_0(i\\omega)` in imaginary frequencies. """)

c.add_property(name = "G_tau",
               getter = cfunction("block_gf_view<imtime> G_tau ()"),
               doc = """Accumulated :math:`G(\\tau)` in imaginary time. """)

c.add_property(name = "G_l",
               getter = cfunction("block_gf_view<legendre> G_l ()"),
               doc = """Accumulated :math:`G_l` in Legendre polynomials representation. """)

c.add_property(name = "G2_iw_inu_inup_pp",
               getter = cfunction("block2_gf_view<cartesian_product<imfreq,imfreq,imfreq>,tensor_valued<4>> G2_iw_inu_inup_pp ()"),
               doc = """Accumulated two-particle Green\xe2\x80\x99s function :math:`G^{(2)}(i\\omega,i\\nu,i\\nu\')` in the pp-channel. """)

c.add_property(name = "G2_iw_inu_inup_ph",
               getter = cfunction("block2_gf_view<cartesian_product<imfreq,imfreq,imfreq>,tensor_valued<4>> G2_iw_inu_inup_ph ()"),
               doc = """Accumulated two-particle Green\xe2\x80\x99s function :math:`G^{(2)}(i\\omega,i\\nu,i\\nu\')` in the ph-channel. """)

c.add_property(name = "G2_iw_l_lp_pp",
               getter = cfunction("block2_gf_view<cartesian_product<imfreq,legendre,legendre>,tensor_valued<4>> G2_iw_l_lp_pp ()"),
               doc = """Accumulated two-particle Green\xe2\x80\x99s function :math:`G^{(2)}(i\\omega,l,l\')` in the pp-channel. """)

c.add_property(name = "G2_iw_l_lp_ph",
               getter = cfunction("block2_gf_view<cartesian_product<imfreq,legendre,legendre>,tensor_valued<4>> G2_iw_l_lp_ph ()"),
               doc = """Accumulated two-particle Green\xe2\x80\x99s function :math:`G^{(2)}(i\\omega,l,l\')` in the ph-channel. """)

c.add_property(name = "atomic_gf",
               getter = cfunction("block_gf_view<imtime> atomic_gf ()"),
               doc = """Atomic :math:`G(\\tau)` in imaginary time. """)

c.add_property(name = "density_matrix",
               getter = cfunction("std::vector<matrix_t> density_matrix ()"),
               doc = """Accumulated density matrix. """)

c.add_property(name = "h_loc_diagonalization",
               getter = cfunction("cthyb::atom_diag h_loc_diagonalization ()"),
               doc = """Diagonalization of :math:`H_{loc}`. """)

c.add_property(name = "perturbation_order_total",
               getter = cfunction("triqs::statistics::histogram get_perturbation_order_total ()"),
               doc = """Histogram of the total perturbation order. """)

c.add_property(name = "perturbation_order",
               getter = cfunction("histo_map_t get_perturbation_order ()"),
               doc = """Histograms of the perturbation order for each block. """)

c.add_property(name = "performance_analysis",
               getter = cfunction("histo_map_t get_performance_analysis ()"),
               doc = """Histograms related to the performance analysis. """)

c.add_property(name = "average_sign",
               getter = cfunction("mc_weight_t average_sign ()"),
               doc = """Monte Carlo average sign. """)

c.add_property(name = "solve_status",
               getter = cfunction("int solve_status ()"),
               doc = """Status of the ``solve()`` on exit. """)

module.add_class(c)

# The class atom_diag
c = class_(
        py_type = "AtomDiag",  # name of the python class
        c_type = "atom_diag",   # name of the C++ class
        hdf5 = True,
        doc = r"",   # doc of the C++ class
)

c.add_constructor("""(many_body_op_t h_, triqs::hilbert_space::fundamental_operator_set fops)""",
                  doc = """ """)

c.add_constructor("""(many_body_op_t h_, triqs::hilbert_space::fundamental_operator_set fops, std::vector<many_body_op_t> qn_vector)""",
                  doc = """ """)

c.add_method("""int get_block_dim (int b)""",
             doc = """The dimension of block b """)

c.add_method("""int flatten_block_index (int block_index, int i)""",
             doc = """Returns the index in the full hilbert space for block_index and i, the index within the block. """)

c.add_method("""double get_eigenvalue (int block_index, int i)""",
             doc = """Get the i-th eigenvalue of block bl """)

c.add_method("""long c_connection (int op_linear_index, int block_index)""",
             doc = """Connections for fundamental operators C\n\n op_linear_index : the linear index (i.e. number) of the c operator, as defined by the fundamental_operator_set fops\n block_number : the number of the initial block\n @return : the number of the final block """)

c.add_method("""long cdag_connection (int op_linear_index, int block_index)""",
             doc = """Connections for fundamental operators C^\\dagger\n\n op_linear_index : the linear index (i.e. number) of the c operator, as defined by the fundamental_operator_set fops\n block_number : the number of the initial block\n @return : the number of the final block """)

c.add_method("""matrix<h_scalar_t> c_matrix (int op_linear_index, int block_index)""",
             doc = """Matrix for fundamental operators C\n\n op_linear_index : the linear index (i.e. number) of the c operator, as defined by the fundamental_operator_set fops\n block_number : the number of the initial block\n @return : the number of the final block """)

c.add_method("""matrix<h_scalar_t> cdag_matrix (int op_linear_index, int block_index)""",
             doc = """Matrix for fundamental operators C^\\dagger\n\n op_linear_index : the linear index (i.e. number) of the c operator, as defined by the fundamental_operator_set fops\n block_number : the number of the initial block\n @return : the number of the final block """)

c.add_property(name = "h_atomic",
               getter = cfunction("many_body_op_t get_h_atomic ()"),
               doc = """The Hamiltonian """)

c.add_property(name = "fops",
               getter = cfunction("triqs::hilbert_space::fundamental_operator_set get_fops ()"),
               doc = """The fundamental operator set used at construction """)

c.add_property(name = "full_hilbert_space_dim",
               getter = cfunction("int get_full_hilbert_space_dim ()"),
               doc = """Dimension of the full Hilbert space """)

c.add_property(name = "n_blocks",
               getter = cfunction("int n_blocks ()"),
               doc = """Number of Blocks """)

c.add_property(name = "fock_states",
               getter = cfunction("std::vector<std::vector<fock_state_t>> get_fock_states ()"),
               doc = """The list of fock states for each block """)

c.add_property(name = "unitary_matrices",
               getter = cfunction("std::vector<matrix<h_scalar_t>> get_unitary_matrices ()"),
               doc = """Unitary matrices that transform from Fock states to atomic eigenstates """)

c.add_property(name = "energies",
               getter = cfunction("std::vector<std::vector<double>> get_energies ()"),
               doc = """A vector of all the energies, by blocks. result[block_number][i] is the energy """)

c.add_property(name = "quantum_numbers",
               getter = cfunction("std::vector<std::vector<double>> get_quantum_numbers ()"),
               doc = """A vector of all the QNs, by blocks : result[block_number][qn_index] is the ..... """)

c.add_property(name = "gs_energy",
               getter = cfunction("double get_gs_energy ()"),
               doc = """Ground state energy (i.e. min of all subspaces). """)

c.add_property(name = "vacuum_block_index",
               getter = cfunction("int get_vacuum_block_index ()"),
               doc = """Returns the block index of the vacuum state. """)

c.add_property(name = "vacuum_inner_index",
               getter = cfunction("int get_vacuum_inner_index ()"),
               doc = """Returns the inner index of the vacuum state. """)

c.add_property(name = "vacuum_state",
               getter = cfunction("full_hilbert_space_state_t get_vacuum_state ()"),
               doc = """Returns the vacuum state as a long vector in the full Hilbert space. """)

module.add_class(c)

#  Free functions

module.add_function ("double partition_function (cthyb::atom_diag atom, double beta)", doc = """The atomic partition function""")

module.add_function ("block_matrix_t atomic_density_matrix (cthyb::atom_diag atom, double beta)", doc = """The atomic density matrix""")

module.add_function ("block_gf<imtime> atomic_gf (cthyb::atom_diag atom, double beta, std::map<std::string,indices_t> indices_list, int n_tau, std::vector<std::pair<int,int>> excluded_states = {})", doc = """The atomic green function, possibly with excluded states (default none)""")

module.add_function ("double trace_rho_op (block_matrix_t density_matrix, many_body_op_t op, cthyb::atom_diag atom)", doc = "Trace (op * density_matrix)")

module.add_function ("full_hilbert_space_state_t act (many_body_op_t op, full_hilbert_space_state_t st, cthyb::atom_diag atom)", doc = """Act with operator op on state st""")

module.add_function ("std::vector<std::vector<double>> quantum_number_eigenvalues (many_body_op_t op, cthyb::atom_diag atom)", doc = """The operator op is supposed to be a quantum number (if not -> exception)\n @return the eigenvalue by block""")

module.add_function ("std::vector<std::vector<double>> quantum_number_eigenvalues2 (many_body_op_t op, cthyb::atom_diag atom)", doc = """The operator op is supposed to be a quantum number (if not -> exception)\n @return the eigenvalue by block""")

module.generate_code()
