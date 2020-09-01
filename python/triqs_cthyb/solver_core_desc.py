# Generated automatically using the command :
# c++2py ../../c++/triqs_cthyb/solver_core.hpp -p --members_read_only -N triqs_cthyb -a triqs_cthyb -m solver_core -o solver_core --moduledoc="The TRIQS cthyb solver" --cxxflags="-std=c++17" -C triqs --only="solver_core block_order"
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "solver_core", doc = r"The TRIQS cthyb solver", app_name = "triqs_cthyb")

# Imports
module.add_imports(*['triqs.atom_diag', 'triqs.gf', 'triqs.operators', 'triqs.statistics.histograms', 'h5._h5py'])

# Add here all includes
module.add_include("triqs_cthyb/solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/map.hpp>
#include <cpp2py/converters/optional.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/set.hpp>
#include <cpp2py/converters/std_array.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/variant.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/arrays.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>

using namespace triqs_cthyb;
""")

module.add_enum("block_order", ['block_order::AABB', 'block_order::ABBA'], "triqs_cthyb", doc = r"""Order of block indices for Block2Gf objects""")

# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "triqs_cthyb::solver_core",   # name of the C++ class
        doc = r"""Core class of the cthyb solver""",   # doc of the C++ class
        hdf5 = True,
)

c.add_member(c_name = "G_tau",
             c_type = "std::optional<G_tau_t>",
             read_only= True,
             doc = r"""Single-particle Green's function :math:`G(\tau)` in imaginary time.""")

c.add_member(c_name = "G_tau_accum",
             c_type = "std::optional<G_tau_G_target_t>",
             read_only= True,
             doc = r"""Intermediate Green's function to accumulate g(tau), either real or complex""")

c.add_member(c_name = "asymmetry_G_tau",
             c_type = "std::optional<G_tau_G_target_t>",
             read_only= True,
             doc = r"""Violation of the fundamental Green function property G(tau)[i,j] = G(tau)*[j,i] after the measurement""")

c.add_member(c_name = "G_l",
             c_type = "std::optional<G_l_t>",
             read_only= True,
             doc = r"""Single-particle Green's function :math:`G_l` in Legendre polynomial representation.""")

c.add_member(c_name = "O_tau",
             c_type = "std::optional<gf<imtime, scalar_valued> >",
             read_only= True,
             doc = r"""General operator Green's function :math:`O(\tau)` in imaginary time.""")

c.add_member(c_name = "G2_tau",
             c_type = "std::optional<G2_tau_t>",
             read_only= True,
             doc = r"""Two-particle Green's function :math:`G^{(2)}(\tau_1,\tau_2,\tau_3)` (three Fermionic imaginary times)""")

c.add_member(c_name = "G2_iw",
             c_type = "std::optional<G2_iw_t>",
             read_only= True,
             doc = r"""Two-particle Green's function :math:`G^{(2)}(i\nu,i\nu',i\nu'')` (three Fermionic frequencies)""")

c.add_member(c_name = "G2_iw_nfft",
             c_type = "std::optional<G2_iw_t>",
             read_only= True,
             doc = r"""Two-particle Green's function :math:`G^{(2)}(i\nu,i\nu',i\nu'')` (three Fermionic frequencies)""")

c.add_member(c_name = "G2_iw_pp",
             c_type = "std::optional<G2_iw_t>",
             read_only= True,
             doc = r"""Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the pp-channel (one bosonic matsubara and two fermionic)""")

c.add_member(c_name = "G2_iw_pp_nfft",
             c_type = "std::optional<G2_iw_t>",
             read_only= True,
             doc = r"""Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the pp-channel (one bosonic matsubara and two fermionic)""")

c.add_member(c_name = "G2_iw_ph",
             c_type = "std::optional<G2_iw_t>",
             read_only= True,
             doc = r"""Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the ph-channel (one bosonic matsubara and two fermionic)""")

c.add_member(c_name = "G2_iw_ph_nfft",
             c_type = "std::optional<G2_iw_t>",
             read_only= True,
             doc = r"""Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the ph-channel (one bosonic matsubara and two fermionic)""")

c.add_member(c_name = "G2_iwll_pp",
             c_type = "std::optional<G2_iwll_t>",
             read_only= True,
             doc = r"""Two-particle Green's function :math:`G^{(2)}(i\omega,l,l')` in the pp-channel (one bosonic matsubara and two legendre)""")

c.add_member(c_name = "G2_iwll_ph",
             c_type = "std::optional<G2_iwll_t>",
             read_only= True,
             doc = r"""Two-particle Green's function :math:`G^{(2)}(i\omega,l,l')` in the ph-channel (one bosonic matsubara and two legendre)""")

c.add_member(c_name = "perturbation_order_total",
             c_type = "std::optional<histogram>",
             read_only= True,
             doc = r"""Histogram of the total perturbation order""")

c.add_member(c_name = "perturbation_order",
             c_type = "std::optional<histo_map_t>",
             read_only= True,
             doc = r"""Histograms of the perturbation order for each block""")

c.add_member(c_name = "constr_parameters",
             c_type = "triqs_cthyb::constr_parameters_t",
             read_only= True,
             doc = r"""""")

c.add_member(c_name = "solve_parameters",
             c_type = "triqs_cthyb::solve_parameters_t",
             read_only= True,
             doc = r"""""")

c.add_constructor("""(**triqs_cthyb::constr_parameters_t)""", doc = r"""Construct a CTHYB solver



+----------------+-----------------------------------+---------+-----------------------------------------------------------------+
| Parameter Name | Type                              | Default | Documentation                                                   |
+================+===================================+=========+=================================================================+
| beta           | double                            | --      | Inverse temperature                                             |
+----------------+-----------------------------------+---------+-----------------------------------------------------------------+
| gf_struct      | triqs::hilbert_space::gf_struct_t | --      | block structure of the gf                                       |
+----------------+-----------------------------------+---------+-----------------------------------------------------------------+
| n_iw           | int                               | 1025    | Number of Matsubara frequencies for gf<imfreq, matrix_valued>   |
+----------------+-----------------------------------+---------+-----------------------------------------------------------------+
| n_tau          | int                               | 10001   | Number of tau points for gf<imtime, matrix_valued>              |
+----------------+-----------------------------------+---------+-----------------------------------------------------------------+
| n_l            | int                               | 50      | Number of Legendre polynomials for gf<legendre, matrix_valued>  |
+----------------+-----------------------------------+---------+-----------------------------------------------------------------+
""")

c.add_method("""void solve (**triqs_cthyb::solve_parameters_t)""",
             doc = r"""Solve method that performs CTHYB calculation



+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| Parameter Name                | Type                                                      | Default                                                   | Documentation                                                                                                     |
+===============================+===========================================================+===========================================================+===================================================================================================================+
| h_int                         | Operator                                                  | --                                                        | Interacting part of the atomic Hamiltonian                                                                        |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| n_cycles                      | int                                                       | --                                                        | Number of QMC cycles                                                                                              |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| partition_method              | str                                                       | "autopartition"                                           | Partition method                                                                                                  |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| quantum_numbers               | list(Operator)                                            | []                                                        | Quantum numbers                                                                                                   |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| loc_n_min                     | int                                                       | 0                                                         | Restrict local Hilbert space to states with at least this number of particles                                     |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| loc_n_max                     | int                                                       | INT_MAX                                                   | Restrict local Hilbert space to states with at most this number of particles                                      |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| length_cycle                  | int                                                       | 50                                                        | Length of a single QMC cycle                                                                                      |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| n_warmup_cycles               | int                                                       | 5000                                                      | Number of cycles for thermalization                                                                               |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| random_seed                   | int                                                       | 34788 + 928374 * MPI.rank                                 | Seed for random number generator                                                                                  |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| random_name                   | str                                                       | ""                                                        | Name of random number generator                                                                                   |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| max_time                      | int                                                       | -1 = infinite                                             | Maximum runtime in seconds, use -1 to set infinite                                                                |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| verbosity                     | int                                                       | 3 on MPI rank 0, 0 otherwise.                             | Verbosity level                                                                                                   |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_shift                    | bool                                                      | true                                                      | Add shifting an operator as a move?                                                                               |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_double                   | bool                                                      | true                                                      | Add double insertions as a move?                                                                                  |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| use_trace_estimator           | bool                                                      | false                                                     | Calculate the full trace or use an estimate?                                                                      |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G_tau                 | bool                                                      | true                                                      | Measure G(tau)? :math:`G_{ij}(\tau)=G_{ji}^*(\tau)` is enforced for the resulting G(tau)                          |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G_l                   | bool                                                      | false                                                     | Measure G_l (Legendre)? Note, no hermiticity in G_l is ensured                                                    |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_O_tau                 | std::optional<std::pair<many_body_op_t, many_body_op_t> > | std::optional<std::pair<many_body_op_t,many_body_op_t>>{} | Measure O_tau by insertion                                                                                        |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_O_tau_min_ins         | int                                                       | 10                                                        | Minumum of operator insertions in: O_tau by insertion measure                                                     |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_tau                | bool                                                      | false                                                     | Measure G^4(tau,tau',tau'') with three fermionic times.                                                           |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_iw                 | bool                                                      | false                                                     | Measure G^4(inu,inu',inu'') with three fermionic frequencies.                                                     |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_iw_nfft            | bool                                                      | false                                                     | Measure G^4(inu,inu',inu'') with three fermionic frequencies.                                                     |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_iw_pp              | bool                                                      | false                                                     | Measure G^4(iomega,inu,inu') within the particle-particle channel.                                                |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_iw_pp_nfft         | bool                                                      | false                                                     | Measure G^4(iomega,inu,inu') within the particle-particle channel.                                                |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_iw_ph              | bool                                                      | false                                                     | Measure G^4(iomega,inu,inu') within the particle-hole channel.                                                    |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_iw_ph_nfft         | bool                                                      | false                                                     | Measure G^4(iomega,inu,inu') within the particle-hole channel.                                                    |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_iwll_pp            | bool                                                      | false                                                     | Measure G^2(iomega,l,l') within the particle-particle channel.                                                    |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_iwll_ph            | bool                                                      | false                                                     | Measure G^2(iomega,l,l') within the particle-hole channel.                                                        |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_block_order        | triqs_cthyb::block_order                                  | block_order::AABB                                         | Order of block indices in the definition of G^2.                                                                  |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_blocks             | std::set<std::pair<std::string, std::string> >            | measure all blocks                                        | List of block index pairs of G^2 to measure.                                                                      |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_n_tau              | int                                                       | 10                                                        | Number of imaginary time slices for G^4 measurement.                                                              |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_n_bosonic          | int                                                       | 30                                                        | Number of bosonic Matsubara frequencies for G^4 measurement.                                                      |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_n_fermionic        | int                                                       | 30                                                        | Number of fermionic Matsubara frequencies for G^4 measurement.                                                    |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_n_l                | int                                                       | 20                                                        | Number of Legendre coefficients for G^4(iomega,l,l') measurement.                                                 |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_G2_iwll_nfft_buf_size | int                                                       | 100                                                       | NFFT buffer size for G^4(iomega,l,l') measurement.                                                                |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| nfft_buf_sizes                | std::map<std::string, int>                                | 100 for every block                                       | NFFT buffer sizes for different blocks                                                                            |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_pert_order            | bool                                                      | false                                                     | Measure perturbation order?                                                                                       |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| measure_density_matrix        | bool                                                      | false                                                     | Measure the reduced impurity density matrix?                                                                      |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| use_norm_as_weight            | bool                                                      | false                                                     | Use the norm of the density matrix in the weight if true, otherwise use Trace                                     |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| performance_analysis          | bool                                                      | false                                                     | Analyse performance of trace computation with histograms (developers only)?                                       |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| proposal_prob                 | dict(str:float)                                           | {}                                                        | Operator insertion/removal probabilities for different blocks                                                     |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_global                   | dict(str : dict(indices : indices))                       | {}                                                        | List of global moves (with their names). Each move is specified with an index substitution dictionary.            |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| move_global_prob              | double                                                    | 0.05                                                      | Overall probability of the global moves                                                                           |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| imag_threshold                | double                                                    | 1.e-13                                                    | Threshold below which imaginary components of Delta and h_loc are set to zero                                     |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| det_init_size                 | int                                                       | 100                                                       | The maximum size of the determinant matrix before a resize                                                        |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| det_n_operations_before_check | int                                                       | 100                                                       | Max number of ops before the test of deviation of the det, M^-1 is performed.                                     |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| det_precision_warning         | double                                                    | 1.e-8                                                     | Threshold for determinant precision warnings                                                                      |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| det_precision_error           | double                                                    | 1.e-5                                                     | Threshold for determinant precision error                                                                         |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
| det_singular_threshold        | double                                                    | -1                                                        | Bound for the determinant matrix being singular, abs(det) > singular_threshold. If <0, it is !isnormal(abs(det))  |
+-------------------------------+-----------------------------------------------------------+-----------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------+
""")

c.add_property(name = "h_loc",
               getter = cfunction("triqs_cthyb::many_body_op_t h_loc ()"),
               doc = r"""The local Hamiltonian of the problem: :math:`H_{loc}` used in the last call to ``solve()``.""")

c.add_property(name = "last_constr_parameters",
               getter = cfunction("triqs_cthyb::constr_parameters_t last_constr_parameters ()"),
               doc = r"""Set of parameters used in the construction of the ``solver_core`` class.""")

c.add_property(name = "last_solve_parameters",
               getter = cfunction("triqs_cthyb::solve_parameters_t last_solve_parameters ()"),
               doc = r"""Set of parameters used in the last call to ``solve()``.""")

c.add_property(name = "Delta_infty",
               getter = cfunction("std::vector<matrix<dcomplex> > Delta_infty ()"),
               doc = r""":math:`G_0^{-1}(i\omega_n = \infty)` in Matsubara Frequency.""")

c.add_property(name = "Delta_tau",
               getter = cfunction("block_gf_view<triqs::gfs::imtime> Delta_tau ()"),
               doc = r""":math:`\Delta(\tau)` in imaginary time.""")

c.add_property(name = "G0_iw",
               getter = cfunction("block_gf_view<triqs::gfs::imfreq> G0_iw ()"),
               doc = r""":math:`G_0(i\omega)` in imaginary frequencies.""")

c.add_property(name = "density_matrix",
               getter = cfunction("std::vector<matrix_t> density_matrix ()"),
               doc = r"""Accumulated density matrix.""")

c.add_property(name = "h_loc_diagonalization",
               getter = cfunction("triqs_cthyb::atom_diag h_loc_diagonalization ()"),
               doc = r"""Diagonalization of :math:`H_{loc}`.""")

c.add_property(name = "performance_analysis",
               getter = cfunction("triqs_cthyb::histo_map_t get_performance_analysis ()"),
               doc = r"""Histograms related to the performance analysis.""")

c.add_property(name = "average_sign",
               getter = cfunction("triqs_cthyb::mc_weight_t average_sign ()"),
               doc = r"""Monte Carlo average sign.""")

c.add_property(name = "solve_status",
               getter = cfunction("int solve_status ()"),
               doc = r"""Status of the ``solve()`` on exit.""")

module.add_class(c)


# Converter for solve_parameters_t
c = converter_(
        c_type = "triqs_cthyb::solve_parameters_t",
        doc = r"""""",
)
c.add_member(c_name = "h_int",
             c_type = "triqs_cthyb::many_body_op_t",
             initializer = """  """,
             doc = r"""Interacting part of the atomic Hamiltonian
     type: Operator""")

c.add_member(c_name = "n_cycles",
             c_type = "int",
             initializer = """  """,
             doc = r"""Number of QMC cycles""")

c.add_member(c_name = "partition_method",
             c_type = "std::string",
             initializer = """ "autopartition" """,
             doc = r"""Partition method
     type: str""")

c.add_member(c_name = "quantum_numbers",
             c_type = "std::vector<many_body_op_t>",
             initializer = """ std::vector<many_body_op_t>{} """,
             doc = r"""Quantum numbers
     type: list(Operator)
     default: []""")

c.add_member(c_name = "loc_n_min",
             c_type = "int",
             initializer = """ 0 """,
             doc = r"""Restrict local Hilbert space to states with at least this number of particles
     default: 0""")

c.add_member(c_name = "loc_n_max",
             c_type = "int",
             initializer = """ INT_MAX """,
             doc = r"""Restrict local Hilbert space to states with at most this number of particles
     default: INT_MAX""")

c.add_member(c_name = "length_cycle",
             c_type = "int",
             initializer = """ 50 """,
             doc = r"""Length of a single QMC cycle
     default: 50""")

c.add_member(c_name = "n_warmup_cycles",
             c_type = "int",
             initializer = """ 5000 """,
             doc = r"""Number of cycles for thermalization
     default: 5000""")

c.add_member(c_name = "random_seed",
             c_type = "int",
             initializer = """ 34788+928374*mpi::communicator().rank() """,
             doc = r"""Seed for random number generator
     default: 34788 + 928374 * MPI.rank""")

c.add_member(c_name = "random_name",
             c_type = "std::string",
             initializer = """ "" """,
             doc = r"""Name of random number generator
     type: str""")

c.add_member(c_name = "max_time",
             c_type = "int",
             initializer = """ -1 """,
             doc = r"""Maximum runtime in seconds, use -1 to set infinite
     default: -1 = infinite""")

c.add_member(c_name = "verbosity",
             c_type = "int",
             initializer = """ ((mpi::communicator().rank()==0)?3:0) """,
             doc = r"""Verbosity level
     default: 3 on MPI rank 0, 0 otherwise.""")

c.add_member(c_name = "move_shift",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Add shifting an operator as a move?""")

c.add_member(c_name = "move_double",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Add double insertions as a move?""")

c.add_member(c_name = "use_trace_estimator",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Calculate the full trace or use an estimate?""")

c.add_member(c_name = "measure_G_tau",
             c_type = "bool",
             initializer = """ true """,
             doc = r"""Measure G(tau)? :math:`G_{ij}(\tau)=G_{ji}^*(\tau)` is enforced for the resulting G(tau)""")

c.add_member(c_name = "measure_G_l",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure G_l (Legendre)? Note, no hermiticity in G_l is ensured""")

c.add_member(c_name = "measure_O_tau",
             c_type = "std::optional<std::pair<many_body_op_t, many_body_op_t> >",
             initializer = """ std::optional<std::pair<many_body_op_t,many_body_op_t>>{} """,
             doc = r"""Measure O_tau by insertion""")

c.add_member(c_name = "measure_O_tau_min_ins",
             c_type = "int",
             initializer = """ 10 """,
             doc = r"""Minumum of operator insertions in: O_tau by insertion measure""")

c.add_member(c_name = "measure_G2_tau",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure G^4(tau,tau',tau'') with three fermionic times.""")

c.add_member(c_name = "measure_G2_iw",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure G^4(inu,inu',inu'') with three fermionic frequencies.""")

c.add_member(c_name = "measure_G2_iw_nfft",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure G^4(inu,inu',inu'') with three fermionic frequencies.""")

c.add_member(c_name = "measure_G2_iw_pp",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure G^4(iomega,inu,inu') within the particle-particle channel.""")

c.add_member(c_name = "measure_G2_iw_pp_nfft",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure G^4(iomega,inu,inu') within the particle-particle channel.""")

c.add_member(c_name = "measure_G2_iw_ph",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure G^4(iomega,inu,inu') within the particle-hole channel.""")

c.add_member(c_name = "measure_G2_iw_ph_nfft",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure G^4(iomega,inu,inu') within the particle-hole channel.""")

c.add_member(c_name = "measure_G2_iwll_pp",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure G^2(iomega,l,l') within the particle-particle channel.""")

c.add_member(c_name = "measure_G2_iwll_ph",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure G^2(iomega,l,l') within the particle-hole channel.""")

c.add_member(c_name = "measure_G2_block_order",
             c_type = "triqs_cthyb::block_order",
             initializer = """ block_order::AABB """,
             doc = r"""Order of block indices in the definition of G^2.""")

c.add_member(c_name = "measure_G2_blocks",
             c_type = "std::set<std::pair<std::string, std::string> >",
             initializer = """ (std::set<std::pair<std::string,std::string>>{}) """,
             doc = r"""List of block index pairs of G^2 to measure.
     default: measure all blocks""")

c.add_member(c_name = "measure_G2_n_tau",
             c_type = "int",
             initializer = """ 10 """,
             doc = r"""Number of imaginary time slices for G^4 measurement.""")

c.add_member(c_name = "measure_G2_n_bosonic",
             c_type = "int",
             initializer = """ 30 """,
             doc = r"""Number of bosonic Matsubara frequencies for G^4 measurement.""")

c.add_member(c_name = "measure_G2_n_fermionic",
             c_type = "int",
             initializer = """ 30 """,
             doc = r"""Number of fermionic Matsubara frequencies for G^4 measurement.""")

c.add_member(c_name = "measure_G2_n_l",
             c_type = "int",
             initializer = """ 20 """,
             doc = r"""Number of Legendre coefficients for G^4(iomega,l,l') measurement.""")

c.add_member(c_name = "measure_G2_iwll_nfft_buf_size",
             c_type = "int",
             initializer = """ 100 """,
             doc = r"""NFFT buffer size for G^4(iomega,l,l') measurement.""")

c.add_member(c_name = "nfft_buf_sizes",
             c_type = "std::map<std::string, int>",
             initializer = """ (std::map<std::string,int>{}) """,
             doc = r"""NFFT buffer sizes for different blocks
     default: 100 for every block""")

c.add_member(c_name = "measure_pert_order",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure perturbation order?""")

c.add_member(c_name = "measure_density_matrix",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Measure the reduced impurity density matrix?""")

c.add_member(c_name = "use_norm_as_weight",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Use the norm of the density matrix in the weight if true, otherwise use Trace""")

c.add_member(c_name = "performance_analysis",
             c_type = "bool",
             initializer = """ false """,
             doc = r"""Analyse performance of trace computation with histograms (developers only)?""")

c.add_member(c_name = "proposal_prob",
             c_type = "std::map<std::string, double>",
             initializer = """ (std::map<std::string,double>{}) """,
             doc = r"""Operator insertion/removal probabilities for different blocks
     type: dict(str:float)
     default: {}""")

c.add_member(c_name = "move_global",
             c_type = "std::map<std::string, indices_map_t>",
             initializer = """ (std::map<std::string,indices_map_t>{}) """,
             doc = r"""List of global moves (with their names).
     Each move is specified with an index substitution dictionary.
     type: dict(str : dict(indices : indices))
     default: {}""")

c.add_member(c_name = "move_global_prob",
             c_type = "double",
             initializer = """ 0.05 """,
             doc = r"""Overall probability of the global moves""")

c.add_member(c_name = "imag_threshold",
             c_type = "double",
             initializer = """ 1.e-13 """,
             doc = r"""Threshold below which imaginary components of Delta and h_loc are set to zero""")

c.add_member(c_name = "det_init_size",
             c_type = "int",
             initializer = """ 100 """,
             doc = r"""The maximum size of the determinant matrix before a resize""")

c.add_member(c_name = "det_n_operations_before_check",
             c_type = "int",
             initializer = """ 100 """,
             doc = r"""Max number of ops before the test of deviation of the det, M^-1 is performed.""")

c.add_member(c_name = "det_precision_warning",
             c_type = "double",
             initializer = """ 1.e-8 """,
             doc = r"""Threshold for determinant precision warnings""")

c.add_member(c_name = "det_precision_error",
             c_type = "double",
             initializer = """ 1.e-5 """,
             doc = r"""Threshold for determinant precision error""")

c.add_member(c_name = "det_singular_threshold",
             c_type = "double",
             initializer = """ -1 """,
             doc = r"""Bound for the determinant matrix being singular, abs(det) > singular_threshold. If <0, it is !isnormal(abs(det))""")

module.add_converter(c)

# Converter for constr_parameters_t
c = converter_(
        c_type = "triqs_cthyb::constr_parameters_t",
        doc = r"""""",
)
c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = r"""Inverse temperature""")

c.add_member(c_name = "gf_struct",
             c_type = "triqs::hilbert_space::gf_struct_t",
             initializer = """  """,
             doc = r"""block structure of the gf""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 1025 """,
             doc = r"""Number of Matsubara frequencies for gf<imfreq, matrix_valued>""")

c.add_member(c_name = "n_tau",
             c_type = "int",
             initializer = """ 10001 """,
             doc = r"""Number of tau points for gf<imtime, matrix_valued>""")

c.add_member(c_name = "n_l",
             c_type = "int",
             initializer = """ 50 """,
             doc = r"""Number of Legendre polynomials for gf<legendre, matrix_valued>""")

module.add_converter(c)


module.generate_code()