
// C.f. https://numpy.org/doc/1.21/reference/c-api/array.html#importing-the-api
#define PY_ARRAY_UNIQUE_SYMBOL _cpp2py_ARRAY_API
#ifndef CLAIR_C2PY_WRAP_GEN
#ifdef __clang__
// #pragma clang diagnostic ignored "-W#warnings"
#endif
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#pragma GCC diagnostic ignored "-Wcpp"
#endif

#define C2PY_VERSION_MAJOR 0
#define C2PY_VERSION_MINOR 1

#include <c2py/c2py.hpp>
#include <c2py/serialization/h5.hpp>

using c2py::operator""_a;

// ==================== enums =====================

template <> constexpr bool c2py::is_wrapped<triqs_cthyb::block_order> = true;
template <>
const std::map<triqs_cthyb::block_order, str_t> c2py::enum_to_string<triqs_cthyb::block_order> = {{triqs_cthyb::block_order::AABB, "AABB"},
                                                                                                  {triqs_cthyb::block_order::ABBA, "ABBA"}};

// ==================== module classes =====================

// --------- class _c2py_cls_0 -----------
using _c2py_cls_0                                            = triqs_cthyb::constr_parameters_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_0>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_0> = "triqs_cthyb.solver_core.ConstrParametersT";

static int synth_constructor_0(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing triqs_cthyb::constr_parameters_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_0> *)self)->_c = new _c2py_cls_0{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError, ("Error in constructing triqs_cthyb::constr_parameters_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  de("beta", self_c.beta, false);
  de("gf_struct", self_c.gf_struct, false);
  de("n_iw", self_c.n_iw, true);
  de("n_tau", self_c.n_tau, true);
  de("n_l", self_c.n_l, true);
  de("delta_interface", self_c.delta_interface, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_0> = synth_constructor_0;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_0> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
beta : {par_0}

gf_struct : {par_1}

n_iw : {par_2}, default=1025

n_tau : {par_3}, default=10001

n_l : {par_4}, default=50

delta_interface : {par_5}, default=false

)DOC",
                      "par",
                      {c2py::python_typename<double>(), c2py::python_typename<triqs::gfs::gf_struct_t>(), c2py::python_typename<int>(),
                       c2py::python_typename<int>(), c2py::python_typename<int>(), c2py::python_typename<bool>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_0>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_0 = R"DOC(Inverse temperature :math:`\beta`.)DOC";
constexpr auto _c2py_doc_member_1 = R"DOC(Structure of the Green's function (names and sizes of blocks).)DOC";
constexpr auto _c2py_doc_member_2 = R"DOC(Number of Matsubara frequencies.)DOC";
constexpr auto _c2py_doc_member_3 = R"DOC(Number of imaginary-time points.)DOC";
constexpr auto _c2py_doc_member_4 = R"DOC(Number of Legendre polynomials.)DOC";
constexpr auto _c2py_doc_member_5 = R"DOC(Use :math:`\Delta(\tau)` and :math:`h_{loc0}` as input instead of :math:`G_0(i\omega)`.)DOC";
static PyObject *prop_get_dict_0(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  c2py::pydict dic;
  dic["beta"]            = self_c.beta;
  dic["gf_struct"]       = self_c.gf_struct;
  dic["n_iw"]            = self_c.n_iw;
  dic["n_tau"]           = self_c.n_tau;
  dic["n_l"]             = self_c.n_l;
  dic["delta_interface"] = self_c.delta_interface;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_0>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_0::beta, _c2py_cls_0>("beta", _c2py_doc_member_0),
   c2py::getsetdef_from_member<&_c2py_cls_0::gf_struct, _c2py_cls_0>("gf_struct", _c2py_doc_member_1),
   c2py::getsetdef_from_member<&_c2py_cls_0::n_iw, _c2py_cls_0>("n_iw", _c2py_doc_member_2),
   c2py::getsetdef_from_member<&_c2py_cls_0::n_tau, _c2py_cls_0>("n_tau", _c2py_doc_member_3),
   c2py::getsetdef_from_member<&_c2py_cls_0::n_l, _c2py_cls_0>("n_l", _c2py_doc_member_4),
   c2py::getsetdef_from_member<&_c2py_cls_0::delta_interface, _c2py_cls_0>("delta_interface", _c2py_doc_member_5),
   {"__dict__", (getter)prop_get_dict_0, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_0> =
   R"DOC(Parameters used for constructing the solver class.)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_0>;
// --------- class _c2py_cls_1 -----------
using _c2py_cls_1                                            = triqs_cthyb::solve_parameters_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_1>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_1> = "triqs_cthyb.solver_core.SolveParametersT";

static int synth_constructor_1(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing triqs_cthyb::solve_parameters_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_1> *)self)->_c = new _c2py_cls_1{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError, ("Error in constructing triqs_cthyb::solve_parameters_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_1> *)self)->_c);
  de("h_int", self_c.h_int, false);
  de("n_cycles", self_c.n_cycles, false);
  de("partition_method", self_c.partition_method, true);
  de("quantum_numbers", self_c.quantum_numbers, true);
  de("loc_n_min", self_c.loc_n_min, true);
  de("loc_n_max", self_c.loc_n_max, true);
  de("length_cycle", self_c.length_cycle, true);
  de("n_warmup_cycles", self_c.n_warmup_cycles, true);
  de("random_seed", self_c.random_seed, true);
  de("random_name", self_c.random_name, true);
  de("max_time", self_c.max_time, true);
  de("verbosity", self_c.verbosity, true);
  de("move_shift", self_c.move_shift, true);
  de("move_double", self_c.move_double, true);
  de("use_trace_estimator", self_c.use_trace_estimator, true);
  de("measure_G_tau", self_c.measure_G_tau, true);
  de("measure_G_l", self_c.measure_G_l, true);
  de("measure_O_tau", self_c.measure_O_tau, true);
  de("measure_F_tau", self_c.measure_F_tau, true);
  de("worm_eta", self_c.worm_eta, true);
  de("worm_prob", self_c.worm_prob, true);
  de("measure_O_tau_min_ins", self_c.measure_O_tau_min_ins, true);
  de("measure_G2_tau", self_c.measure_G2_tau, true);
  de("measure_G2_iw", self_c.measure_G2_iw, true);
  de("measure_G2_iw_nfft", self_c.measure_G2_iw_nfft, true);
  de("measure_G2_iw_pp", self_c.measure_G2_iw_pp, true);
  de("measure_G2_iw_pp_nfft", self_c.measure_G2_iw_pp_nfft, true);
  de("measure_G2_iw_ph", self_c.measure_G2_iw_ph, true);
  de("measure_G2_iw_ph_nfft", self_c.measure_G2_iw_ph_nfft, true);
  de("measure_G2_iwll_pp", self_c.measure_G2_iwll_pp, true);
  de("measure_G2_iwll_ph", self_c.measure_G2_iwll_ph, true);
  de("measure_G2_block_order", self_c.measure_G2_block_order, true);
  de("measure_G2_blocks", self_c.measure_G2_blocks, true);
  de("measure_G2_n_tau", self_c.measure_G2_n_tau, true);
  de("measure_G2_n_bosonic", self_c.measure_G2_n_bosonic, true);
  de("measure_G2_n_fermionic", self_c.measure_G2_n_fermionic, true);
  de("measure_G2_n_l", self_c.measure_G2_n_l, true);
  de("measure_G2_iwll_nfft_buf_size", self_c.measure_G2_iwll_nfft_buf_size, true);
  de("nfft_buf_sizes", self_c.nfft_buf_sizes, true);
  de("measure_pert_order", self_c.measure_pert_order, true);
  de("measure_density_matrix", self_c.measure_density_matrix, true);
  de("use_norm_as_weight", self_c.use_norm_as_weight, true);
  de("initial_configuration", self_c.initial_configuration, true);
  de("performance_analysis", self_c.performance_analysis, true);
  de("proposal_prob", self_c.proposal_prob, true);
  de("move_global", self_c.move_global, true);
  de("move_global_prob", self_c.move_global_prob, true);
  de("imag_threshold", self_c.imag_threshold, true);
  de("det_init_size", self_c.det_init_size, true);
  de("det_n_operations_before_check", self_c.det_n_operations_before_check, true);
  de("det_precision_warning", self_c.det_precision_warning, true);
  de("det_precision_error", self_c.det_precision_error, true);
  de("det_singular_threshold", self_c.det_singular_threshold, true);
  de("off_diag_threshold", self_c.off_diag_threshold, true);
  de("h_loc0", self_c.h_loc0, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_1> = synth_constructor_1;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_1> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
h_int : {par_0}

n_cycles : {par_1}

partition_method : {par_2}, default="autopartition"

quantum_numbers : {par_3}, default={}

loc_n_min : {par_4}, default=0

loc_n_max : {par_5}, default=INT_MAX

length_cycle : {par_6}, default=50

n_warmup_cycles : {par_7}, default=5000

random_seed : {par_8}, default=34788 + 928374 * mpi::communicator().rank()

random_name : {par_9}, default=""

max_time : {par_10}, default=-1

verbosity : {par_11}, default== 0) ? 3 : 0)

move_shift : {par_12}, default=true

move_double : {par_13}, default=true

use_trace_estimator : {par_14}, default=false

measure_G_tau : {par_15}, default=true

measure_G_l : {par_16}, default=false

measure_O_tau : {par_17}, default={}

measure_O_tau_min_ins : {par_18}, default=10

measure_G2_tau : {par_19}, default=false

measure_G2_iw : {par_20}, default=false

measure_G2_iw_nfft : {par_21}, default=false

measure_G2_iw_pp : {par_22}, default=false

measure_G2_iw_pp_nfft : {par_23}, default=false

measure_G2_iw_ph : {par_24}, default=false

measure_G2_iw_ph_nfft : {par_25}, default=false

measure_G2_iwll_pp : {par_26}, default=false

measure_G2_iwll_ph : {par_27}, default=false

measure_G2_block_order : {par_28}, default=block_order::AABB

measure_G2_blocks : {par_29}, default={}

measure_G2_n_tau : {par_30}, default=10

measure_G2_n_bosonic : {par_31}, default=30

measure_G2_n_fermionic : {par_32}, default=30

measure_G2_n_l : {par_33}, default=20

measure_G2_iwll_nfft_buf_size : {par_34}, default=100

nfft_buf_sizes : {par_35}, default={}

measure_pert_order : {par_36}, default=false

measure_density_matrix : {par_37}, default=false

use_norm_as_weight : {par_38}, default=false

initial_configuration : {par_39}, default={}

performance_analysis : {par_40}, default=false

proposal_prob : {par_41}, default={}

move_global : {par_42}, default={}

move_global_prob : {par_43}, default=0.05

imag_threshold : {par_44}, default=1.e-13

det_init_size : {par_45}, default=100

det_n_operations_before_check : {par_46}, default=100

det_precision_warning : {par_47}, default=1.e-8

det_precision_error : {par_48}, default=1.e-5

det_singular_threshold : {par_49}, default=-1

off_diag_threshold : {par_50}, default=0.0

h_loc0 : {par_51}, default={}

measure_F_tau : {par_52}, default=false

worm_eta : {par_53}, default=1.0

worm_prob : {par_54}, default=0.3

)DOC",
                      "par",
                      {c2py::python_typename<triqs_cthyb::many_body_op_t>(),
                       c2py::python_typename<long>(),
                       c2py::python_typename<std::string>(),
                       c2py::python_typename<std::vector<triqs_cthyb::many_body_op_t>>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<long>(),
                       c2py::python_typename<long>(),
                       c2py::python_typename<long>(),
                       c2py::python_typename<std::string>(),
                       c2py::python_typename<long>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<std::optional<std::pair<triqs_cthyb::many_body_op_t, triqs_cthyb::many_body_op_t>>>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<triqs_cthyb::block_order>(),
                       c2py::python_typename<std::set<std::pair<std::string, std::string>>>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<std::map<std::string, long>>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<std::optional<triqs_cthyb::configuration>>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<std::map<std::string, double>>(),
                       c2py::python_typename<std::map<std::string, triqs_cthyb::indices_map_t>>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<int>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<std::optional<triqs_cthyb::many_body_op_t>>(),
                       c2py::python_typename<bool>(),
                       c2py::python_typename<double>(),
                       c2py::python_typename<double>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_1>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_6  = R"DOC(Interacting part of the atomic Hamiltonian.)DOC";
constexpr auto _c2py_doc_member_7  = R"DOC(Number of QMC cycles.)DOC";
constexpr auto _c2py_doc_member_8  = R"DOC(Partition method.)DOC";
constexpr auto _c2py_doc_member_9  = R"DOC(Quantum numbers.)DOC";
constexpr auto _c2py_doc_member_10 = R"DOC(Restrict local Hilbert space to states with at least this number of particles.)DOC";
constexpr auto _c2py_doc_member_11 = R"DOC(Restrict local Hilbert space to states with at most this number of particles.)DOC";
constexpr auto _c2py_doc_member_12 = R"DOC(Length of a single QMC cycle.)DOC";
constexpr auto _c2py_doc_member_13 = R"DOC(Number of cycles for thermalization.)DOC";
constexpr auto _c2py_doc_member_14 = R"DOC(Seed for random number generator.)DOC";
constexpr auto _c2py_doc_member_15 = R"DOC(Name of random number generator.)DOC";
constexpr auto _c2py_doc_member_16 = R"DOC(Maximum runtime in seconds, use -1 to set infinite.)DOC";
constexpr auto _c2py_doc_member_17 = R"DOC(Verbosity level.)DOC";
constexpr auto _c2py_doc_member_18 = R"DOC(Add shifting an operator as a move?)DOC";
constexpr auto _c2py_doc_member_19 = R"DOC(Add double insertions as a move?)DOC";
constexpr auto _c2py_doc_member_20 = R"DOC(Calculate the full trace or use an estimate?)DOC";
constexpr auto _c2py_doc_member_21 = R"DOC(Measure :math:`G(\tau)`? Hermiticity :math:`G_{ij}(\tau) = G_{ji}^*(\tau)` is enforced.)DOC";
constexpr auto _c2py_doc_member_22 = R"DOC(Measure :math:`G_l` (Legendre)? No hermiticity is enforced.)DOC";
constexpr auto _c2py_doc_member_23 = R"DOC(Measure :math:`O(\tau)` by insertion.)DOC";
constexpr auto _c2py_doc_member_76 = R"DOC(Measure the improved-estimator correlator :math:`F(\tau)` by worm sampling.)DOC";
constexpr auto _c2py_doc_member_77 = R"DOC(Extended-ensemble weight for the :math:`F(\tau)` worm sector.)DOC";
constexpr auto _c2py_doc_member_78 = R"DOC(Relative proposal weight for :math:`F(\tau)` worm insert/remove/shift moves.)DOC";
constexpr auto _c2py_doc_member_24 = R"DOC(Minimum number of operator insertions in the :math:`O(\tau)` insertion measure.)DOC";
constexpr auto _c2py_doc_member_25 = R"DOC(Measure :math:`G^{(2)}(\tau,\tau',\tau'')` with three fermionic times.)DOC";
constexpr auto _c2py_doc_member_26 = R"DOC(Measure :math:`G^{(2)}(i\nu,i\nu',i\nu'')` with three fermionic frequencies.)DOC";
constexpr auto _c2py_doc_member_27 = R"DOC(Measure :math:`G^{(2)}(i\nu,i\nu',i\nu'')` with three fermionic frequencies.)DOC";
constexpr auto _c2py_doc_member_28 = R"DOC(Measure :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-particle channel.)DOC";
constexpr auto _c2py_doc_member_29 = R"DOC(Measure :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-particle channel.)DOC";
constexpr auto _c2py_doc_member_30 = R"DOC(Measure :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-hole channel.)DOC";
constexpr auto _c2py_doc_member_31 = R"DOC(Measure :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-hole channel.)DOC";
constexpr auto _c2py_doc_member_32 = R"DOC(Measure :math:`G^{(2)}(i\omega,l,l')` in the particle-particle channel.)DOC";
constexpr auto _c2py_doc_member_33 = R"DOC(Measure :math:`G^{(2)}(i\omega,l,l')` in the particle-hole channel.)DOC";
constexpr auto _c2py_doc_member_34 = R"DOC(Order of block indices in the definition of :math:`G^{(2)}`.)DOC";
constexpr auto _c2py_doc_member_35 = R"DOC(List of block index pairs of :math:`G^{(2)}` to measure.)DOC";
constexpr auto _c2py_doc_member_36 = R"DOC(Number of imaginary-time slices for the :math:`G^{(2)}` measurement.)DOC";
constexpr auto _c2py_doc_member_37 = R"DOC(Number of bosonic Matsubara frequencies for the :math:`G^{(2)}` measurement.)DOC";
constexpr auto _c2py_doc_member_38 = R"DOC(Number of fermionic Matsubara frequencies for the :math:`G^{(2)}` measurement.)DOC";
constexpr auto _c2py_doc_member_39 = R"DOC(Number of Legendre coefficients for the :math:`G^{(2)}(i\omega,l,l')` measurement.)DOC";
constexpr auto _c2py_doc_member_40 = R"DOC(NFFT buffer size for the :math:`G^{(2)}(i\omega,l,l')` measurement.)DOC";
constexpr auto _c2py_doc_member_41 = R"DOC(NFFT buffer sizes for different blocks.)DOC";
constexpr auto _c2py_doc_member_42 = R"DOC(Measure perturbation order?)DOC";
constexpr auto _c2py_doc_member_43 = R"DOC(Measure the reduced impurity density matrix?)DOC";
constexpr auto _c2py_doc_member_44 = R"DOC(Use the norm of the density matrix in the weight (instead of the trace)?)DOC";
constexpr auto _c2py_doc_member_45 = R"DOC(Initial configuration of the run (advanced, use with care).)DOC";
constexpr auto _c2py_doc_member_46 = R"DOC(Analyse performance of the trace computation with histograms (developers only)?)DOC";
constexpr auto _c2py_doc_member_47 = R"DOC(Operator insertion/removal probabilities for different blocks.)DOC";
constexpr auto _c2py_doc_member_48 =
   R"DOC(List of global moves (with their names). Each move is specified with an index substitution dictionary.)DOC";
constexpr auto _c2py_doc_member_49 = R"DOC(Overall probability of the global moves.)DOC";
constexpr auto _c2py_doc_member_50 = R"DOC(Threshold below which imaginary components of :math:`\Delta` and :math:`h_{loc}` are set to zero.)DOC";
constexpr auto _c2py_doc_member_51 = R"DOC(The maximum size of the determinant matrix before a resize.)DOC";
constexpr auto _c2py_doc_member_52 = R"DOC(Maximum number of operations before testing the accuracy of :math:`\det(M)` and :math:`M^{-1}`.)DOC";
constexpr auto _c2py_doc_member_53 = R"DOC(Threshold for determinant precision warnings.)DOC";
constexpr auto _c2py_doc_member_54 = R"DOC(Threshold for determinant precision error.)DOC";
constexpr auto _c2py_doc_member_55 = R"DOC(Bound for the determinant matrix being singular (if :math:`< 0`, checks for subnormal numbers).)DOC";
constexpr auto _c2py_doc_member_56 = R"DOC(Threshold below which off-diagonal components of :math:`h_{loc}` are set to zero.)DOC";
constexpr auto _c2py_doc_member_57 = R"DOC(Quadratic part of the local Hamiltonian. Must be provided if the :math:`\Delta` interface is used.)DOC";
static PyObject *prop_get_dict_1(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_1> *)self)->_c);
  c2py::pydict dic;
  dic["h_int"]                         = self_c.h_int;
  dic["n_cycles"]                      = self_c.n_cycles;
  dic["partition_method"]              = self_c.partition_method;
  dic["quantum_numbers"]               = self_c.quantum_numbers;
  dic["loc_n_min"]                     = self_c.loc_n_min;
  dic["loc_n_max"]                     = self_c.loc_n_max;
  dic["length_cycle"]                  = self_c.length_cycle;
  dic["n_warmup_cycles"]               = self_c.n_warmup_cycles;
  dic["random_seed"]                   = self_c.random_seed;
  dic["random_name"]                   = self_c.random_name;
  dic["max_time"]                      = self_c.max_time;
  dic["verbosity"]                     = self_c.verbosity;
  dic["move_shift"]                    = self_c.move_shift;
  dic["move_double"]                   = self_c.move_double;
  dic["use_trace_estimator"]           = self_c.use_trace_estimator;
  dic["measure_G_tau"]                 = self_c.measure_G_tau;
  dic["measure_G_l"]                   = self_c.measure_G_l;
  dic["measure_O_tau"]                 = self_c.measure_O_tau;
  dic["measure_F_tau"]                 = self_c.measure_F_tau;
  dic["worm_eta"]                      = self_c.worm_eta;
  dic["worm_prob"]                     = self_c.worm_prob;
  dic["measure_O_tau_min_ins"]         = self_c.measure_O_tau_min_ins;
  dic["measure_G2_tau"]                = self_c.measure_G2_tau;
  dic["measure_G2_iw"]                 = self_c.measure_G2_iw;
  dic["measure_G2_iw_nfft"]            = self_c.measure_G2_iw_nfft;
  dic["measure_G2_iw_pp"]              = self_c.measure_G2_iw_pp;
  dic["measure_G2_iw_pp_nfft"]         = self_c.measure_G2_iw_pp_nfft;
  dic["measure_G2_iw_ph"]              = self_c.measure_G2_iw_ph;
  dic["measure_G2_iw_ph_nfft"]         = self_c.measure_G2_iw_ph_nfft;
  dic["measure_G2_iwll_pp"]            = self_c.measure_G2_iwll_pp;
  dic["measure_G2_iwll_ph"]            = self_c.measure_G2_iwll_ph;
  dic["measure_G2_block_order"]        = self_c.measure_G2_block_order;
  dic["measure_G2_blocks"]             = self_c.measure_G2_blocks;
  dic["measure_G2_n_tau"]              = self_c.measure_G2_n_tau;
  dic["measure_G2_n_bosonic"]          = self_c.measure_G2_n_bosonic;
  dic["measure_G2_n_fermionic"]        = self_c.measure_G2_n_fermionic;
  dic["measure_G2_n_l"]                = self_c.measure_G2_n_l;
  dic["measure_G2_iwll_nfft_buf_size"] = self_c.measure_G2_iwll_nfft_buf_size;
  dic["nfft_buf_sizes"]                = self_c.nfft_buf_sizes;
  dic["measure_pert_order"]            = self_c.measure_pert_order;
  dic["measure_density_matrix"]        = self_c.measure_density_matrix;
  dic["use_norm_as_weight"]            = self_c.use_norm_as_weight;
  dic["initial_configuration"]         = self_c.initial_configuration;
  dic["performance_analysis"]          = self_c.performance_analysis;
  dic["proposal_prob"]                 = self_c.proposal_prob;
  dic["move_global"]                   = self_c.move_global;
  dic["move_global_prob"]              = self_c.move_global_prob;
  dic["imag_threshold"]                = self_c.imag_threshold;
  dic["det_init_size"]                 = self_c.det_init_size;
  dic["det_n_operations_before_check"] = self_c.det_n_operations_before_check;
  dic["det_precision_warning"]         = self_c.det_precision_warning;
  dic["det_precision_error"]           = self_c.det_precision_error;
  dic["det_singular_threshold"]        = self_c.det_singular_threshold;
  dic["off_diag_threshold"]            = self_c.off_diag_threshold;
  dic["h_loc0"]                        = self_c.h_loc0;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_1>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_1::h_int, _c2py_cls_1>("h_int", _c2py_doc_member_6),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_cycles, _c2py_cls_1>("n_cycles", _c2py_doc_member_7),
   c2py::getsetdef_from_member<&_c2py_cls_1::partition_method, _c2py_cls_1>("partition_method", _c2py_doc_member_8),
   c2py::getsetdef_from_member<&_c2py_cls_1::quantum_numbers, _c2py_cls_1>("quantum_numbers", _c2py_doc_member_9),
   c2py::getsetdef_from_member<&_c2py_cls_1::loc_n_min, _c2py_cls_1>("loc_n_min", _c2py_doc_member_10),
   c2py::getsetdef_from_member<&_c2py_cls_1::loc_n_max, _c2py_cls_1>("loc_n_max", _c2py_doc_member_11),
   c2py::getsetdef_from_member<&_c2py_cls_1::length_cycle, _c2py_cls_1>("length_cycle", _c2py_doc_member_12),
   c2py::getsetdef_from_member<&_c2py_cls_1::n_warmup_cycles, _c2py_cls_1>("n_warmup_cycles", _c2py_doc_member_13),
   c2py::getsetdef_from_member<&_c2py_cls_1::random_seed, _c2py_cls_1>("random_seed", _c2py_doc_member_14),
   c2py::getsetdef_from_member<&_c2py_cls_1::random_name, _c2py_cls_1>("random_name", _c2py_doc_member_15),
   c2py::getsetdef_from_member<&_c2py_cls_1::max_time, _c2py_cls_1>("max_time", _c2py_doc_member_16),
   c2py::getsetdef_from_member<&_c2py_cls_1::verbosity, _c2py_cls_1>("verbosity", _c2py_doc_member_17),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_shift, _c2py_cls_1>("move_shift", _c2py_doc_member_18),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_double, _c2py_cls_1>("move_double", _c2py_doc_member_19),
   c2py::getsetdef_from_member<&_c2py_cls_1::use_trace_estimator, _c2py_cls_1>("use_trace_estimator", _c2py_doc_member_20),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G_tau, _c2py_cls_1>("measure_G_tau", _c2py_doc_member_21),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G_l, _c2py_cls_1>("measure_G_l", _c2py_doc_member_22),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_O_tau, _c2py_cls_1>("measure_O_tau", _c2py_doc_member_23),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_F_tau, _c2py_cls_1>("measure_F_tau", _c2py_doc_member_76),
   c2py::getsetdef_from_member<&_c2py_cls_1::worm_eta, _c2py_cls_1>("worm_eta", _c2py_doc_member_77),
   c2py::getsetdef_from_member<&_c2py_cls_1::worm_prob, _c2py_cls_1>("worm_prob", _c2py_doc_member_78),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_O_tau_min_ins, _c2py_cls_1>("measure_O_tau_min_ins", _c2py_doc_member_24),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_tau, _c2py_cls_1>("measure_G2_tau", _c2py_doc_member_25),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_iw, _c2py_cls_1>("measure_G2_iw", _c2py_doc_member_26),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_iw_nfft, _c2py_cls_1>("measure_G2_iw_nfft", _c2py_doc_member_27),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_iw_pp, _c2py_cls_1>("measure_G2_iw_pp", _c2py_doc_member_28),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_iw_pp_nfft, _c2py_cls_1>("measure_G2_iw_pp_nfft", _c2py_doc_member_29),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_iw_ph, _c2py_cls_1>("measure_G2_iw_ph", _c2py_doc_member_30),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_iw_ph_nfft, _c2py_cls_1>("measure_G2_iw_ph_nfft", _c2py_doc_member_31),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_iwll_pp, _c2py_cls_1>("measure_G2_iwll_pp", _c2py_doc_member_32),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_iwll_ph, _c2py_cls_1>("measure_G2_iwll_ph", _c2py_doc_member_33),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_block_order, _c2py_cls_1>("measure_G2_block_order", _c2py_doc_member_34),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_blocks, _c2py_cls_1>("measure_G2_blocks", _c2py_doc_member_35),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_n_tau, _c2py_cls_1>("measure_G2_n_tau", _c2py_doc_member_36),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_n_bosonic, _c2py_cls_1>("measure_G2_n_bosonic", _c2py_doc_member_37),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_n_fermionic, _c2py_cls_1>("measure_G2_n_fermionic", _c2py_doc_member_38),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_n_l, _c2py_cls_1>("measure_G2_n_l", _c2py_doc_member_39),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_G2_iwll_nfft_buf_size, _c2py_cls_1>("measure_G2_iwll_nfft_buf_size", _c2py_doc_member_40),
   c2py::getsetdef_from_member<&_c2py_cls_1::nfft_buf_sizes, _c2py_cls_1>("nfft_buf_sizes", _c2py_doc_member_41),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_pert_order, _c2py_cls_1>("measure_pert_order", _c2py_doc_member_42),
   c2py::getsetdef_from_member<&_c2py_cls_1::measure_density_matrix, _c2py_cls_1>("measure_density_matrix", _c2py_doc_member_43),
   c2py::getsetdef_from_member<&_c2py_cls_1::use_norm_as_weight, _c2py_cls_1>("use_norm_as_weight", _c2py_doc_member_44),
   c2py::getsetdef_from_member<&_c2py_cls_1::initial_configuration, _c2py_cls_1>("initial_configuration", _c2py_doc_member_45),
   c2py::getsetdef_from_member<&_c2py_cls_1::performance_analysis, _c2py_cls_1>("performance_analysis", _c2py_doc_member_46),
   c2py::getsetdef_from_member<&_c2py_cls_1::proposal_prob, _c2py_cls_1>("proposal_prob", _c2py_doc_member_47),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_global, _c2py_cls_1>("move_global", _c2py_doc_member_48),
   c2py::getsetdef_from_member<&_c2py_cls_1::move_global_prob, _c2py_cls_1>("move_global_prob", _c2py_doc_member_49),
   c2py::getsetdef_from_member<&_c2py_cls_1::imag_threshold, _c2py_cls_1>("imag_threshold", _c2py_doc_member_50),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_init_size, _c2py_cls_1>("det_init_size", _c2py_doc_member_51),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_n_operations_before_check, _c2py_cls_1>("det_n_operations_before_check", _c2py_doc_member_52),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_precision_warning, _c2py_cls_1>("det_precision_warning", _c2py_doc_member_53),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_precision_error, _c2py_cls_1>("det_precision_error", _c2py_doc_member_54),
   c2py::getsetdef_from_member<&_c2py_cls_1::det_singular_threshold, _c2py_cls_1>("det_singular_threshold", _c2py_doc_member_55),
   c2py::getsetdef_from_member<&_c2py_cls_1::off_diag_threshold, _c2py_cls_1>("off_diag_threshold", _c2py_doc_member_56),
   c2py::getsetdef_from_member<&_c2py_cls_1::h_loc0, _c2py_cls_1>("h_loc0", _c2py_doc_member_57),
   {"__dict__", (getter)prop_get_dict_1, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_1> =
   R"DOC(Parameters passed to the solve method of the solver class.)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_1>;
// --------- class _c2py_cls_2 -----------
using _c2py_cls_2                                            = triqs_cthyb::solver_core;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_2>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_2> = "triqs_cthyb.solver_core.SolverCore";
static const auto _c2py_init_0 = c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_2, const triqs_cthyb::constr_parameters_t &>("p")};
template <> constexpr initproc c2py::tp_init<_c2py_cls_2> = c2py::pyfkw_constructor<_c2py_init_0>;
template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_2> = _c2py_init_0.doc(R"DOC(
Construct a CTHYB solver.

Parameters
----------
p : {par_0}
   Parameters used for constructing the solver.
)DOC",
                                                                    {{c2py::python_typename<const triqs_cthyb::constr_parameters_t &>()}});
// solve
static auto const _c2py_fun_0 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_2 &self, const triqs_cthyb::solve_parameters_t &p) -> decltype(auto) { return self.solve(p); }, "self", "p")};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(R"DOC(
Solve the impurity problem.

Parameters
----------
p : {par_0}
   Parameters controlling the Monte Carlo simulation and measurements.
)DOC",
                                                {{c2py::python_typename<const triqs_cthyb::solve_parameters_t &>()}});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_2>[] = {
   {"solve", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_2>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_2>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_2>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_58 = R"DOC(Parameters used for constructing the solver.)DOC";
constexpr auto _c2py_doc_member_59 = R"DOC(Parameters passed to the solve method.)DOC";
constexpr auto _c2py_doc_member_60 = R"DOC(Single-particle Green's function :math:`G(\tau)` in imaginary time.)DOC";
constexpr auto _c2py_doc_member_61 = R"DOC(Intermediate Green's function used to accumulate :math:`G(\tau)` (real or complex).)DOC";
constexpr auto _c2py_doc_member_62 = R"DOC(Violation of the property :math:`G_{ij}(\tau) = G_{ji}^*(\tau)` after the measurement.)DOC";
constexpr auto _c2py_doc_member_63 = R"DOC(Single-particle Green's function :math:`G_l` in the Legendre representation.)DOC";
constexpr auto _c2py_doc_member_64 = R"DOC(General operator Green's function :math:`O(\tau)` in imaginary time.)DOC";
constexpr auto _c2py_doc_member_79 = R"DOC(Improved-estimator correlator :math:`F(\tau)` in imaginary time.)DOC";
constexpr auto _c2py_doc_member_80 = R"DOC(Intermediate Green's function used to accumulate :math:`F(\tau)` (real or complex).)DOC";
constexpr auto _c2py_doc_member_65 = R"DOC(Two-particle Green's function :math:`G^{(2)}(\tau_1,\tau_2,\tau_3)` with three fermionic times.)DOC";
constexpr auto _c2py_doc_member_66 = R"DOC(Two-particle Green's function :math:`G^{(2)}(i\nu,i\nu',i\nu'')` with three fermionic frequencies.)DOC";
constexpr auto _c2py_doc_member_67 = R"DOC(Two-particle Green's function :math:`G^{(2)}(i\nu,i\nu',i\nu'')` with three fermionic frequencies.)DOC";
constexpr auto _c2py_doc_member_68 = R"DOC(Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-particle channel.)DOC";
constexpr auto _c2py_doc_member_69 = R"DOC(Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-particle channel.)DOC";
constexpr auto _c2py_doc_member_70 = R"DOC(Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-hole channel.)DOC";
constexpr auto _c2py_doc_member_71 = R"DOC(Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-hole channel.)DOC";
constexpr auto _c2py_doc_member_72 = R"DOC(Two-particle Green's function :math:`G^{(2)}(i\omega,l,l')` in the particle-particle channel.)DOC";
constexpr auto _c2py_doc_member_73 = R"DOC(Two-particle Green's function :math:`G^{(2)}(i\omega,l,l')` in the particle-hole channel.)DOC";
constexpr auto _c2py_doc_member_74 = R"DOC(Histogram of the total perturbation order.)DOC";
constexpr auto _c2py_doc_member_75 = R"DOC(Histograms of the perturbation order for each block.)DOC";
static constexpr auto prop_doc_0   = R"DOC(:math:`G_0^{-1}(i\omega_n = \infty)` in Matsubara frequencies.)DOC";
static constexpr auto prop_doc_1   = R"DOC(Hybridization function :math:`\Delta(\tau)` in imaginary time.)DOC";
static constexpr auto prop_doc_2   = R"DOC(Non-interacting Green's function :math:`G_0(i\omega)` in Matsubara frequencies.)DOC";
static constexpr auto prop_doc_3   = R"DOC(Auto-correlation time in units of MC cycles.)DOC";
static constexpr auto prop_doc_4 = R"DOC(Whether the auto-correlation time estimate has saturated (false: it is only a lower bound, run longer).)DOC";
static constexpr auto prop_doc_5 = R"DOC(Average perturbation order.)DOC";
static constexpr auto prop_doc_6 = R"DOC(Monte Carlo average sign.)DOC";
static constexpr auto prop_doc_7 = R"DOC(Accumulated density matrix.)DOC";
static constexpr auto prop_doc_8 = R"DOC(The local Hamiltonian :math:`H_{loc}` used in the last solve.)DOC";
static constexpr auto prop_doc_9 = R"DOC(The noninteracting part of the local Hamiltonian.)DOC";
static constexpr auto prop_doc_10 = R"DOC(Diagonalization of :math:`H_{loc}`.)DOC";
static constexpr auto prop_doc_11 = R"DOC(Is the solver compiled with support for complex hybridization?)DOC";
static constexpr auto prop_doc_12 = R"DOC(Final configuration of the last solve call.)DOC";
static constexpr auto prop_doc_13 = R"DOC(Parameters used for constructing the solver.)DOC";
static constexpr auto prop_doc_14 = R"DOC(Parameters used in the last solve.)DOC";
static constexpr auto prop_doc_15 = R"DOC(Is the solver compiled with support for a complex local Hamiltonian?)DOC";
static constexpr auto prop_doc_16 = R"DOC(Histograms related to the performance analysis.)DOC";
static constexpr auto prop_doc_17 = R"DOC(Status of the solve on exit.)DOC";

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_2>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_2::constr_parameters, _c2py_cls_2>("constr_parameters", _c2py_doc_member_58),
   c2py::getsetdef_from_member<&_c2py_cls_2::solve_parameters, _c2py_cls_2>("solve_parameters", _c2py_doc_member_59),
   c2py::getsetdef_from_member<&_c2py_cls_2::G_tau, _c2py_cls_2>("G_tau", _c2py_doc_member_60),
   c2py::getsetdef_from_member<&_c2py_cls_2::G_tau_accum, _c2py_cls_2>("G_tau_accum", _c2py_doc_member_61),
   c2py::getsetdef_from_member<&_c2py_cls_2::asymmetry_G_tau, _c2py_cls_2>("asymmetry_G_tau", _c2py_doc_member_62),
   c2py::getsetdef_from_member<&_c2py_cls_2::G_l, _c2py_cls_2>("G_l", _c2py_doc_member_63),
   c2py::getsetdef_from_member<&_c2py_cls_2::O_tau, _c2py_cls_2>("O_tau", _c2py_doc_member_64),
   c2py::getsetdef_from_member<&_c2py_cls_2::F_tau, _c2py_cls_2>("F_tau", _c2py_doc_member_79),
   c2py::getsetdef_from_member<&_c2py_cls_2::F_tau_accum, _c2py_cls_2>("F_tau_accum", _c2py_doc_member_80),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2_tau, _c2py_cls_2>("G2_tau", _c2py_doc_member_65),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2_iw, _c2py_cls_2>("G2_iw", _c2py_doc_member_66),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2_iw_nfft, _c2py_cls_2>("G2_iw_nfft", _c2py_doc_member_67),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2_iw_pp, _c2py_cls_2>("G2_iw_pp", _c2py_doc_member_68),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2_iw_pp_nfft, _c2py_cls_2>("G2_iw_pp_nfft", _c2py_doc_member_69),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2_iw_ph, _c2py_cls_2>("G2_iw_ph", _c2py_doc_member_70),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2_iw_ph_nfft, _c2py_cls_2>("G2_iw_ph_nfft", _c2py_doc_member_71),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2_iwll_pp, _c2py_cls_2>("G2_iwll_pp", _c2py_doc_member_72),
   c2py::getsetdef_from_member<&_c2py_cls_2::G2_iwll_ph, _c2py_cls_2>("G2_iwll_ph", _c2py_doc_member_73),
   c2py::getsetdef_from_member<&_c2py_cls_2::perturbation_order_total, _c2py_cls_2>("perturbation_order_total", _c2py_doc_member_74),
   c2py::getsetdef_from_member<&_c2py_cls_2::perturbation_order, _c2py_cls_2>("perturbation_order", _c2py_doc_member_75),
   {"Delta_infty", c2py::getter_from_method<c2py::castm<>(&triqs_cthyb::solver_core::Delta_infty)>, nullptr, prop_doc_0, nullptr},
   {"Delta_tau", c2py::getter_from_method<c2py::castm<>(&triqs_cthyb::solver_core::Delta_tau)>, nullptr, prop_doc_1, nullptr},
   {"G0_iw", c2py::getter_from_method<c2py::castm<>(&triqs_cthyb::solver_core::G0_iw)>, nullptr, prop_doc_2, nullptr},
   {"auto_corr_time", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::auto_corr_time)>, nullptr, prop_doc_3, nullptr},
   {"auto_corr_time_converged", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::auto_corr_time_converged)>, nullptr, prop_doc_4,
    nullptr},
   {"average_order", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::average_order)>, nullptr, prop_doc_5, nullptr},
   {"average_sign", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::average_sign)>, nullptr, prop_doc_6, nullptr},
   {"density_matrix", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::density_matrix)>, nullptr, prop_doc_7, nullptr},
   {"h_loc", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::h_loc)>, nullptr, prop_doc_8, nullptr},
   {"h_loc0", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::h_loc0)>, nullptr, prop_doc_9, nullptr},
   {"h_loc_diagonalization", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::h_loc_diagonalization)>, nullptr, prop_doc_10,
    nullptr},
   {"hybridisation_is_complex", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::hybridisation_is_complex)>, nullptr, prop_doc_11,
    nullptr},
   {"last_configuration", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::last_configuration)>, nullptr, prop_doc_12, nullptr},
   {"last_constr_parameters", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::last_constr_parameters)>, nullptr, prop_doc_13,
    nullptr},
   {"last_solve_parameters", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::last_solve_parameters)>, nullptr, prop_doc_14,
    nullptr},
   {"local_hamiltonian_is_complex", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::local_hamiltonian_is_complex)>, nullptr,
    prop_doc_15, nullptr},
   {"performance_analysis", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::get_performance_analysis)>, nullptr, prop_doc_16,
    nullptr},
   {"solve_status", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::solve_status)>, nullptr, prop_doc_17, nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_2> = R"DOC(Continuous-time hybridization-expansion quantum Monte Carlo solver.)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_2>;

// ==================== module functions ====================

//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {PyModuleDef_HEAD_INIT,
                                        "solver_core",                           /* name of module */
                                        R"RAWDOC(The TRIQS cthyb solver)RAWDOC", /* module documentation, may be NULL */
                                        -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
                                        module_methods,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_solver_core() {

  if (not c2py::check_python_version("solver_core")) return NULL;

  // import numpy iff 'numpy/arrayobject.h' included
#ifdef Py_ARRAYOBJECT_H
  import_array();
#endif

  PyObject *m;

  if (PyType_Ready(&c2py::wrap_pytype<c2py::py_range>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_0>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_1>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_2>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  auto &conv_table = *c2py::conv_table_sptr.get();

  conv_table[std::type_index(typeid(c2py::py_range)).name()] = &c2py::wrap_pytype<c2py::py_range>;
#define _add_type(T, N) c2py::add_type_object_to_main<T>(N, m, conv_table)
  _add_type(_c2py_cls_0, "ConstrParametersT");
  _add_type(_c2py_cls_1, "SolveParametersT");
  _add_type(_c2py_cls_2, "SolverCore");
#undef _add_type

  c2py::pyref module = c2py::pyref::module("h5.formats");
  if (not module) return nullptr;
  c2py::pyref register_class = module.attr("register_class");

  register_h5_type<_c2py_cls_2>(register_class);

  return m;
}
#endif
// CLAIR_WRAP_GEN
