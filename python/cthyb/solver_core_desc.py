# Generated automatically using the command :
# c++2py ../../cthyb/solver_core.hpp -p -a cthyb -m solver_core -o solver_core --moduledoc="The cthyb solver" --cxxflags="-std=c++17" -C pytriqs --only="solver_core block_order" -N cthyb -I../../app_atom_diag --include /opt/local/include --include /opt/local/include/openmpi-clang50/
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "solver_core", doc = "The cthyb solver", app_name = "cthyb")

# Imports
import pytriqs.gf
import pytriqs.operators
import pytriqs.statistics.histograms
import pytriqs.atom_diag

# Add here all includes
module.add_include("../../cthyb/solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/map.hpp>
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/set.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/variant.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>

using namespace cthyb;
""")

module.add_enum("block_order", ['block_order::AABB', 'block_order::ABBA'], "cthyb", """Order of block indices for Block2Gf objects""")

# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "cthyb::solver_core",   # name of the C++ class
        doc = """Core class of the cthyb solver""",   # doc of the C++ class
        hdf5 = True,
)

c.add_member(c_name = "G_tau",
             c_type = "std::optional<G_tau_t>",
             read_only= False,
             doc = """Single-particle Green\'s function :math:`G(\\tau)` in imaginary time.""")

c.add_member(c_name = "G_tau_accum",
             c_type = "std::optional<G_tau_G_target_t>",
             read_only= False,
             doc = """Intermediate Green\'s function to accumulate g(tau), either real or complex""")

c.add_member(c_name = "G_l",
             c_type = "std::optional<G_l_t>",
             read_only= False,
             doc = """Single-particle Green\'s function :math:`G_l` in Legendre polynomial representation.""")

c.add_member(c_name = "G2_tau",
             c_type = "std::optional<G2_tau_t>",
             read_only= False,
             doc = """Two-particle Green\'s function :math:`G^{(2)}(\\tau_1,\\tau_2,\\tau_3)` (three Fermionic imaginary times)""")

c.add_member(c_name = "G2_iw",
             c_type = "std::optional<G2_iw_t>",
             read_only= False,
             doc = """Two-particle Green\'s function :math:`G^{(2)}(i\\nu,i\\nu\',i\\nu\'\')` (three Fermionic frequencies)""")

c.add_member(c_name = "G2_iw_pp",
             c_type = "std::optional<G2_iw_t>",
             read_only= False,
             doc = """Two-particle Green\'s function :math:`G^{(2)}(i\\omega,i\\nu,i\\nu\')` in the pp-channel (one bosonic matsubara and two fermionic)""")

c.add_member(c_name = "G2_iw_ph",
             c_type = "std::optional<G2_iw_t>",
             read_only= False,
             doc = """Two-particle Green\'s function :math:`G^{(2)}(i\\omega,i\\nu,i\\nu\')` in the ph-channel (one bosonic matsubara and two fermionic)""")

c.add_member(c_name = "G2_iwll_pp",
             c_type = "std::optional<G2_iwll_t>",
             read_only= False,
             doc = """Two-particle Green\'s function :math:`G^{(2)}(i\\omega,l,l\')` in the pp-channel (one bosonic matsubara and two legendre)""")

c.add_member(c_name = "G2_iwll_ph",
             c_type = "std::optional<G2_iwll_t>",
             read_only= False,
             doc = """Two-particle Green\'s function :math:`G^{(2)}(i\\omega,l,l\')` in the ph-channel (one bosonic matsubara and two legendre)""")

c.add_member(c_name = "constr_parameters",
             c_type = "cthyb::constr_parameters_t",
             read_only= False,
             doc = """""")

c.add_member(c_name = "solve_parameters",
             c_type = "cthyb::solve_parameters_t",
             read_only= False,
             doc = """""")

c.add_constructor("""(**cthyb::constr_parameters_t)""", doc = """Construct a CTHYB solver\n\n :param p: Set of parameters specific to the CTHYB solver
+----------------+--------------------+---------+-----------------------------------------------------------------+
| Parameter Name | Type               | Default | Documentation                                                   |
+================+====================+=========+=================================================================+
| beta           | double             |         | Inverse temperature                                             |
+----------------+--------------------+---------+-----------------------------------------------------------------+
| gf_struct      | cthyb::gf_struct_t |         | block structure of the gf                                       |
+----------------+--------------------+---------+-----------------------------------------------------------------+
| n_iw           | int                | 1025    | Number of Matsubara frequencies for gf<imfreq, matrix_valued>   |
+----------------+--------------------+---------+-----------------------------------------------------------------+
| n_tau          | int                | 10001   | Number of tau points for gf<imtime, matrix_valued>              |
+----------------+--------------------+---------+-----------------------------------------------------------------+
| n_l            | int                | 50      | Number of Legendre polynomials for gf<legendre, matrix_valued>  |
+----------------+--------------------+---------+-----------------------------------------------------------------+""")

c.add_method("""void solve (**cthyb::solve_parameters_t)""",
             doc = """Solve method that performs CTHYB calculation\n\n :param p: Set of parameters for the CTHYB calculation
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Parameter Name                | Type                                           | Default                                          | Documentation                                                                                                                                                                   |
+===============================+================================================+==================================================+=================================================================================================================================================================================+
| h_int                         | cthyb::many_body_op_t                          |                                                  | Interacting part of the atomic Hamiltonian\n     type: Operator                                                                                                                 |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| n_cycles                      | int                                            |                                                  | Number of QMC cycles                                                                                                                                                            |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| partition_method              | std::string                                    | "autopartition"                                  | Partition method\n     type: str                                                                                                                                                |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| quantum_numbers               | std::vector<many_body_op_t>                    | std::vector<many_body_op_t>{}                    | Quantum numbers\n     type: list(Operator)\n     default: []                                                                                                                    |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| length_cycle                  | int                                            | 50                                               | Length of a single QMC cycle\n     default: 50                                                                                                                                  |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| n_warmup_cycles               | int                                            | 5000                                             | Number of cycles for thermalization\n     default: 5000                                                                                                                         |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| random_seed                   | int                                            | 34788+928374*triqs::mpi::communicator().rank()   | Seed for random number generator\n     default: 34788 + 928374 * MPI.rank                                                                                                       |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| random_name                   | std::string                                    | ""                                               | Name of random number generator\n     type: str                                                                                                                                 |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| max_time                      | int                                            | -1                                               | Maximum runtime in seconds, use -1 to set infinite\n     default: -1 = infinite                                                                                                 |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| verbosity                     | int                                            | ((triqs::mpi::communicator().rank()==0)?3:0)     | Verbosity level\n     default: 3 on MPI rank 0, 0 otherwise.                                                                                                                    |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| move_shift                    | bool                                           | true                                             | Add shifting an operator as a move?                                                                                                                                             |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| move_double                   | bool                                           | false                                            | Add double insertions as a move?                                                                                                                                                |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| use_trace_estimator           | bool                                           | false                                            | Calculate the full trace or use an estimate?                                                                                                                                    |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G_tau                 | bool                                           | true                                             | Measure G(tau)?                                                                                                                                                                 |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G_l                   | bool                                           | false                                            | Measure G_l (Legendre)?                                                                                                                                                         |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G2_tau                | bool                                           | false                                            | Measure G^4(tau,tau\',tau\'\') with three fermionic times.                                                                                                                      |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G2_iw                 | bool                                           | false                                            | Measure G^4(inu,inu\',inu\'\') with three fermionic frequencies.                                                                                                                |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G2_iw_pp              | bool                                           | false                                            | Measure G^4(iomega,inu,inu\') within the particle-particle channel.                                                                                                             |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G2_iw_ph              | bool                                           | false                                            | Measure G^4(iomega,inu,inu\') within the particle-hole channel.                                                                                                                 |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G2_iwll_pp            | bool                                           | false                                            | Measure G^2(iomega,l,l\') within the particle-particle channel.                                                                                                                 |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G2_iwll_ph            | bool                                           | false                                            | Measure G^2(iomega,l,l\') within the particle-hole channel.                                                                                                                     |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G2_block_order        | cthyb::block_order                             | block_order::AABB                                | Order of block indices in the definition of G^2.                                                                                                                                |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G2_blocks             | std::set<std::pair<std::string, std::string> > | (std::set<std::pair<std::string,std::string>>{}) | List of block index pairs of G^2 to measure.\n     default: measure all blocks                                                                                                  |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G2_n_tau              | int                                            | 10                                               | Number of imaginary time slices for G^4 measurement.                                                                                                                            |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G2_n_bosonic          | int                                            | 30                                               | Number of bosonic Matsubara frequencies for G^4 measurement.                                                                                                                    |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G2_n_fermionic        | int                                            | 30                                               | Number of fermionic Matsubara frequencies for G^4 measurement.                                                                                                                  |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G2_n_l                | int                                            | 20                                               | Number of Legendre coefficients for G^4(iomega,l,l\') measurement.                                                                                                              |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_G2_iwll_nfft_buf_size | int                                            | 100                                              | NFFT buffer size for G^4(iomega,l,l\') measurement.                                                                                                                             |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| nfft_buf_sizes                | std::map<std::string, int>                     | (std::map<std::string,int>{})                    | NFFT buffer sizes for different blocks\n     default: 100 for every block                                                                                                       |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_pert_order            | bool                                           | false                                            | Measure perturbation order?                                                                                                                                                     |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| measure_density_matrix        | bool                                           | false                                            | Measure the reduced impurity density matrix?                                                                                                                                    |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| use_norm_as_weight            | bool                                           | false                                            | Use the norm of the density matrix in the weight if true, otherwise use Trace                                                                                                   |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| performance_analysis          | bool                                           | false                                            | Analyse performance of trace computation with histograms (developers only)?                                                                                                     |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| proposal_prob                 | std::map<std::string, double>                  | (std::map<std::string,double>{})                 | Operator insertion/removal probabilities for different blocks\n     type: dict(str:float)\n     default: {}                                                                     |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| move_global                   | std::map<std::string, indices_map_t>           | (std::map<std::string,indices_map_t>{})          | List of global moves (with their names).\n     Each move is specified with an index substitution dictionary.\n     type: dict(str : dict(indices : indices))\n     default: {}  |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| move_global_prob              | double                                         | 0.05                                             | Overall probability of the global moves                                                                                                                                         |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| imag_threshold                | double                                         | 1.e-15                                           | Threshold below which imaginary components of Delta and h_loc are set to zero                                                                                                   |
+-------------------------------+------------------------------------------------+--------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+""")

c.add_method("""std::string hdf5_scheme ()""",
             is_static = True,
             doc = """""")

c.add_method("""cthyb::solver_core h5_read_construct (triqs::h5::group h5group, std::string subgroup_name)""",
             is_static = True,
             doc = """""")

c.add_property(name = "h_loc",
               getter = cfunction("cthyb::many_body_op_t h_loc ()"),
               doc = """The local Hamiltonian of the problem: :math:`H_{loc}` used in the last call to ``solve()``.""")

c.add_property(name = "last_constr_parameters",
               getter = cfunction("cthyb::constr_parameters_t last_constr_parameters ()"),
               doc = """Set of parameters used in the construction of the ``solver_core`` class.""")

c.add_property(name = "last_solve_parameters",
               getter = cfunction("cthyb::solve_parameters_t last_solve_parameters ()"),
               doc = """Set of parameters used in the last call to ``solve()``.""")

c.add_property(name = "Delta_tau",
               getter = cfunction("block_gf_view<triqs::gfs::imtime> Delta_tau ()"),
               doc = """:math:`\\Delta(\\tau)` in imaginary time.""")

c.add_property(name = "G0_iw",
               getter = cfunction("block_gf_view<triqs::gfs::imfreq> G0_iw ()"),
               doc = """:math:`G_0(i\\omega)` in imaginary frequencies.""")

c.add_property(name = "density_matrix",
               getter = cfunction("std::vector<matrix_t> density_matrix ()"),
               doc = """Accumulated density matrix.""")

c.add_property(name = "h_loc_diagonalization",
               getter = cfunction("cthyb::atom_diag h_loc_diagonalization ()"),
               doc = """Diagonalization of :math:`H_{loc}`.""")

c.add_property(name = "perturbation_order_total",
               getter = cfunction("triqs::statistics::histogram get_perturbation_order_total ()"),
               doc = """Histogram of the total perturbation order.""")

c.add_property(name = "perturbation_order",
               getter = cfunction("cthyb::histo_map_t get_perturbation_order ()"),
               doc = """Histograms of the perturbation order for each block.""")

c.add_property(name = "performance_analysis",
               getter = cfunction("cthyb::histo_map_t get_performance_analysis ()"),
               doc = """Histograms related to the performance analysis.""")

c.add_property(name = "average_sign",
               getter = cfunction("cthyb::mc_weight_t average_sign ()"),
               doc = """Monte Carlo average sign.""")

c.add_property(name = "solve_status",
               getter = cfunction("int solve_status ()"),
               doc = """Status of the ``solve()`` on exit.""")

module.add_class(c)


# Converter for solve_parameters_t
c = converter_(
        c_type = "cthyb::solve_parameters_t",
        doc = """""",
)
c.add_member(c_name = "h_int",
             c_type = "cthyb::many_body_op_t",
             initializer = """  """,
             doc = """Interacting part of the atomic Hamiltonian\n     type: Operator""")

c.add_member(c_name = "n_cycles",
             c_type = "int",
             initializer = """  """,
             doc = """Number of QMC cycles""")

c.add_member(c_name = "partition_method",
             c_type = "std::string",
             initializer = """ "autopartition" """,
             doc = """Partition method\n     type: str""")

c.add_member(c_name = "quantum_numbers",
             c_type = "std::vector<many_body_op_t>",
             initializer = """ std::vector<many_body_op_t>{} """,
             doc = """Quantum numbers\n     type: list(Operator)\n     default: []""")

c.add_member(c_name = "length_cycle",
             c_type = "int",
             initializer = """ 50 """,
             doc = """Length of a single QMC cycle\n     default: 50""")

c.add_member(c_name = "n_warmup_cycles",
             c_type = "int",
             initializer = """ 5000 """,
             doc = """Number of cycles for thermalization\n     default: 5000""")

c.add_member(c_name = "random_seed",
             c_type = "int",
             initializer = """ 34788+928374*triqs::mpi::communicator().rank() """,
             doc = """Seed for random number generator\n     default: 34788 + 928374 * MPI.rank""")

c.add_member(c_name = "random_name",
             c_type = "std::string",
             initializer = """ "" """,
             doc = """Name of random number generator\n     type: str""")

c.add_member(c_name = "max_time",
             c_type = "int",
             initializer = """ -1 """,
             doc = """Maximum runtime in seconds, use -1 to set infinite\n     default: -1 = infinite""")

c.add_member(c_name = "verbosity",
             c_type = "int",
             initializer = """ ((triqs::mpi::communicator().rank()==0)?3:0) """,
             doc = """Verbosity level\n     default: 3 on MPI rank 0, 0 otherwise.""")

c.add_member(c_name = "move_shift",
             c_type = "bool",
             initializer = """ true """,
             doc = """Add shifting an operator as a move?""")

c.add_member(c_name = "move_double",
             c_type = "bool",
             initializer = """ false """,
             doc = """Add double insertions as a move?""")

c.add_member(c_name = "use_trace_estimator",
             c_type = "bool",
             initializer = """ false """,
             doc = """Calculate the full trace or use an estimate?""")

c.add_member(c_name = "measure_G_tau",
             c_type = "bool",
             initializer = """ true """,
             doc = """Measure G(tau)?""")

c.add_member(c_name = "measure_G_l",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure G_l (Legendre)?""")

c.add_member(c_name = "measure_G2_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure G^4(tau,tau\',tau\'\') with three fermionic times.""")

c.add_member(c_name = "measure_G2_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure G^4(inu,inu\',inu\'\') with three fermionic frequencies.""")

c.add_member(c_name = "measure_G2_iw_pp",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure G^4(iomega,inu,inu\') within the particle-particle channel.""")

c.add_member(c_name = "measure_G2_iw_ph",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure G^4(iomega,inu,inu\') within the particle-hole channel.""")

c.add_member(c_name = "measure_G2_iwll_pp",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure G^2(iomega,l,l\') within the particle-particle channel.""")

c.add_member(c_name = "measure_G2_iwll_ph",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure G^2(iomega,l,l\') within the particle-hole channel.""")

c.add_member(c_name = "measure_G2_block_order",
             c_type = "cthyb::block_order",
             initializer = """ block_order::AABB """,
             doc = """Order of block indices in the definition of G^2.""")

c.add_member(c_name = "measure_G2_blocks",
             c_type = "std::set<std::pair<std::string, std::string> >",
             initializer = """ (std::set<std::pair<std::string,std::string>>{}) """,
             doc = """List of block index pairs of G^2 to measure.\n     default: measure all blocks""")

c.add_member(c_name = "measure_G2_n_tau",
             c_type = "int",
             initializer = """ 10 """,
             doc = """Number of imaginary time slices for G^4 measurement.""")

c.add_member(c_name = "measure_G2_n_bosonic",
             c_type = "int",
             initializer = """ 30 """,
             doc = """Number of bosonic Matsubara frequencies for G^4 measurement.""")

c.add_member(c_name = "measure_G2_n_fermionic",
             c_type = "int",
             initializer = """ 30 """,
             doc = """Number of fermionic Matsubara frequencies for G^4 measurement.""")

c.add_member(c_name = "measure_G2_n_l",
             c_type = "int",
             initializer = """ 20 """,
             doc = """Number of Legendre coefficients for G^4(iomega,l,l\') measurement.""")

c.add_member(c_name = "measure_G2_iwll_nfft_buf_size",
             c_type = "int",
             initializer = """ 100 """,
             doc = """NFFT buffer size for G^4(iomega,l,l\') measurement.""")

c.add_member(c_name = "nfft_buf_sizes",
             c_type = "std::map<std::string, int>",
             initializer = """ (std::map<std::string,int>{}) """,
             doc = """NFFT buffer sizes for different blocks\n     default: 100 for every block""")

c.add_member(c_name = "measure_pert_order",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure perturbation order?""")

c.add_member(c_name = "measure_density_matrix",
             c_type = "bool",
             initializer = """ false """,
             doc = """Measure the reduced impurity density matrix?""")

c.add_member(c_name = "use_norm_as_weight",
             c_type = "bool",
             initializer = """ false """,
             doc = """Use the norm of the density matrix in the weight if true, otherwise use Trace""")

c.add_member(c_name = "performance_analysis",
             c_type = "bool",
             initializer = """ false """,
             doc = """Analyse performance of trace computation with histograms (developers only)?""")

c.add_member(c_name = "proposal_prob",
             c_type = "std::map<std::string, double>",
             initializer = """ (std::map<std::string,double>{}) """,
             doc = """Operator insertion/removal probabilities for different blocks\n     type: dict(str:float)\n     default: {}""")

c.add_member(c_name = "move_global",
             c_type = "std::map<std::string, indices_map_t>",
             initializer = """ (std::map<std::string,indices_map_t>{}) """,
             doc = """List of global moves (with their names).\n     Each move is specified with an index substitution dictionary.\n     type: dict(str : dict(indices : indices))\n     default: {}""")

c.add_member(c_name = "move_global_prob",
             c_type = "double",
             initializer = """ 0.05 """,
             doc = """Overall probability of the global moves""")

c.add_member(c_name = "imag_threshold",
             c_type = "double",
             initializer = """ 1.e-15 """,
             doc = """Threshold below which imaginary components of Delta and h_loc are set to zero""")

module.add_converter(c)

# Converter for constr_parameters_t
c = converter_(
        c_type = "cthyb::constr_parameters_t",
        doc = """""",
)
c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = """Inverse temperature""")

c.add_member(c_name = "gf_struct",
             c_type = "cthyb::gf_struct_t",
             initializer = """  """,
             doc = """block structure of the gf""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 1025 """,
             doc = """Number of Matsubara frequencies for gf<imfreq, matrix_valued>""")

c.add_member(c_name = "n_tau",
             c_type = "int",
             initializer = """ 10001 """,
             doc = """Number of tau points for gf<imtime, matrix_valued>""")

c.add_member(c_name = "n_l",
             c_type = "int",
             initializer = """ 50 """,
             doc = """Number of Legendre polynomials for gf<legendre, matrix_valued>""")

module.add_converter(c)


module.generate_code()
