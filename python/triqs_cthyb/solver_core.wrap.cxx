
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

#define C2PY_VERSION_MAJOR 1
#define C2PY_VERSION_MINOR 0

#include <c2py/c2py.hpp>
#include <c2py/version_check.hpp>
#include <c2py/serialization/h5.hpp>

using c2py::operator""_a;

// ==================== enums =====================

template <> constexpr bool c2py::is_wrapped<triqs_cthyb::block_order> = true;
template <>
const std::map<triqs_cthyb::block_order, std::string> c2py::enum_to_string<triqs_cthyb::block_order> = {{triqs_cthyb::block_order::AABB, "AABB"},
                                                                                                        {triqs_cthyb::block_order::ABBA, "ABBA"}};

// ==================== module classes =====================

// --------- class _c2py_cls_e57534ae -----------
using _c2py_cls_e57534ae                                            = triqs_cthyb::constr_parameters_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_e57534ae>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_e57534ae> = "triqs_cthyb.solver_core.ConstrParametersT";

static int synth_constructor_99469753(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing triqs_cthyb::constr_parameters_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_e57534ae> *)self)->_c = new _c2py_cls_e57534ae{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError, ("Error in constructing triqs_cthyb::constr_parameters_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_e57534ae> *)self)->_c);
  de("beta", self_c.beta, false);
  de("gf_struct", self_c.gf_struct, false);
  de("n_iw", self_c.n_iw, true);
  de("n_tau", self_c.n_tau, true);
  de("n_l", self_c.n_l, true);
  de("delta_interface", self_c.delta_interface, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_e57534ae> = synth_constructor_99469753;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_e57534ae> =
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
// clang-format off
template <> PyMethodDef c2py::tp_methods<_c2py_cls_e57534ae>[] = {
   {nullptr, nullptr, 0, nullptr} // Sentinel
};
// clang-format on

constexpr auto _c2py_doc_member_b3cf4a40 = R"DOC(Inverse temperature :math:`\beta`.)DOC";
constexpr auto _c2py_doc_member_c3d7f82d = R"DOC(Structure of the Green's function (names and sizes of blocks).)DOC";
constexpr auto _c2py_doc_member_3ad1678d = R"DOC(Number of Matsubara frequencies.)DOC";
constexpr auto _c2py_doc_member_7e4a921d = R"DOC(Number of imaginary-time points.)DOC";
constexpr auto _c2py_doc_member_8c372e07 = R"DOC(Number of Legendre polynomials.)DOC";
constexpr auto _c2py_doc_member_2e016d4c = R"DOC(Use :math:`\Delta(\tau)` and :math:`h_{loc0}` as input instead of :math:`G_0(i\omega)`.)DOC";
static PyObject *prop_get_dict_99469753(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_e57534ae> *)self)->_c);
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
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_e57534ae>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_e57534ae::beta, _c2py_cls_e57534ae>("beta", _c2py_doc_member_b3cf4a40),
   c2py::getsetdef_from_member<&_c2py_cls_e57534ae::gf_struct, _c2py_cls_e57534ae>("gf_struct", _c2py_doc_member_c3d7f82d),
   c2py::getsetdef_from_member<&_c2py_cls_e57534ae::n_iw, _c2py_cls_e57534ae>("n_iw", _c2py_doc_member_3ad1678d),
   c2py::getsetdef_from_member<&_c2py_cls_e57534ae::n_tau, _c2py_cls_e57534ae>("n_tau", _c2py_doc_member_7e4a921d),
   c2py::getsetdef_from_member<&_c2py_cls_e57534ae::n_l, _c2py_cls_e57534ae>("n_l", _c2py_doc_member_8c372e07),
   c2py::getsetdef_from_member<&_c2py_cls_e57534ae::delta_interface, _c2py_cls_e57534ae>("delta_interface", _c2py_doc_member_2e016d4c),
   {"__dict__", (getter)prop_get_dict_99469753, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_e57534ae> =
   R"DOC(Parameters used for constructing the solver class.)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_e57534ae>;
// --------- class _c2py_cls_a2b1dce8 -----------
using _c2py_cls_a2b1dce8                                            = triqs_cthyb::solve_parameters_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_a2b1dce8>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_a2b1dce8> = "triqs_cthyb.solver_core.SolveParametersT";

static int synth_constructor_167a5c32(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing triqs_cthyb::solve_parameters_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_a2b1dce8> *)self)->_c = new _c2py_cls_a2b1dce8{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError, ("Error in constructing triqs_cthyb::solve_parameters_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_a2b1dce8> *)self)->_c);
  de("h_int", self_c.h_int, false);
  de("n_cycles", self_c.n_cycles, false);
  de("partition_method", self_c.partition_method, true);
  de("quantum_numbers", self_c.quantum_numbers, true);
  de("loc_n_min", self_c.loc_n_min, true);
  de("loc_n_max", self_c.loc_n_max, true);
  de("length_cycle", self_c.length_cycle, true);
  de("max_length_cycle", self_c.max_length_cycle, true);
  de("target_auto_corr_time", self_c.target_auto_corr_time, true);
  de("n_warmup_cycles", self_c.n_warmup_cycles, true);
  de("max_warmup_cycles", self_c.max_warmup_cycles, true);
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
  de("measure_densities", self_c.measure_densities, true);
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

template <> constexpr initproc c2py::tp_init<_c2py_cls_a2b1dce8> = synth_constructor_167a5c32;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_a2b1dce8> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
h_int : {par_0}

n_cycles : {par_1}

partition_method : {par_2}, default="autopartition"

quantum_numbers : {par_3}, default={}

loc_n_min : {par_4}, default=0

loc_n_max : {par_5}, default=INT_MAX

length_cycle : {par_6}, default=-1

max_length_cycle : {par_7}, default=5000

target_auto_corr_time : {par_8}, default=2.0

n_warmup_cycles : {par_9}, default=-1

max_warmup_cycles : {par_10}, default=100000

random_seed : {par_11}, default=34788 + 928374 * mpi::communicator().rank()

random_name : {par_12}, default=""

max_time : {par_13}, default=-1

verbosity : {par_14}, default== 0) ? 3 : 0)

move_shift : {par_15}, default=true

move_double : {par_16}, default=true

use_trace_estimator : {par_17}, default=false

measure_G_tau : {par_18}, default=true

measure_G_l : {par_19}, default=false

measure_O_tau : {par_20}, default={}

measure_O_tau_min_ins : {par_21}, default=10

measure_G2_tau : {par_22}, default=false

measure_G2_iw : {par_23}, default=false

measure_G2_iw_nfft : {par_24}, default=false

measure_G2_iw_pp : {par_25}, default=false

measure_G2_iw_pp_nfft : {par_26}, default=false

measure_G2_iw_ph : {par_27}, default=false

measure_G2_iw_ph_nfft : {par_28}, default=false

measure_G2_iwll_pp : {par_29}, default=false

measure_G2_iwll_ph : {par_30}, default=false

measure_G2_block_order : {par_31}, default=block_order::AABB

measure_G2_blocks : {par_32}, default={}

measure_G2_n_tau : {par_33}, default=10

measure_G2_n_bosonic : {par_34}, default=30

measure_G2_n_fermionic : {par_35}, default=30

measure_G2_n_l : {par_36}, default=20

measure_G2_iwll_nfft_buf_size : {par_37}, default=100

nfft_buf_sizes : {par_38}, default={}

measure_densities : {par_39}, default=false

measure_pert_order : {par_40}, default=false

measure_density_matrix : {par_41}, default=false

use_norm_as_weight : {par_42}, default=false

initial_configuration : {par_43}, default={}

performance_analysis : {par_44}, default=false

proposal_prob : {par_45}, default={}

move_global : {par_46}, default={}

move_global_prob : {par_47}, default=0.05

imag_threshold : {par_48}, default=1.e-13

det_init_size : {par_49}, default=100

det_n_operations_before_check : {par_50}, default=100

det_precision_warning : {par_51}, default=1.e-8

det_precision_error : {par_52}, default=1.e-5

det_singular_threshold : {par_53}, default=-1

off_diag_threshold : {par_54}, default=0.0

h_loc0 : {par_55}, default={}

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
                       c2py::python_typename<double>(),
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
                       c2py::python_typename<std::optional<triqs_cthyb::many_body_op_t>>()});

// ----- Method table ----
// clang-format off
template <> PyMethodDef c2py::tp_methods<_c2py_cls_a2b1dce8>[] = {
   {nullptr, nullptr, 0, nullptr} // Sentinel
};
// clang-format on

constexpr auto _c2py_doc_member_19337922 = R"DOC(Interacting part of the atomic Hamiltonian.)DOC";
constexpr auto _c2py_doc_member_e0206678 = R"DOC(Number of QMC cycles.)DOC";
constexpr auto _c2py_doc_member_2a8cd03e = R"DOC(Partition method.)DOC";
constexpr auto _c2py_doc_member_7ae6201c = R"DOC(Quantum numbers.)DOC";
constexpr auto _c2py_doc_member_9a5219d0 = R"DOC(Restrict local Hilbert space to states with at least this number of particles.)DOC";
constexpr auto _c2py_doc_member_903e825a = R"DOC(Restrict local Hilbert space to states with at most this number of particles.)DOC";
constexpr auto _c2py_doc_member_b60dd7e1 = R"DOC(Length of a single QMC cycle (-1: automatically determined from the autocorrelation time).)DOC";
constexpr auto _c2py_doc_member_d9e01446 = R"DOC(Maximum allowed length_cycle when auto-determined (safety cap).)DOC";
constexpr auto _c2py_doc_member_19b5e8ac = R"DOC(Target autocorrelation time in units of length_cycle (used when length_cycle=-1).)DOC";
constexpr auto _c2py_doc_member_261489c3 = R"DOC(Number of cycles for thermalization (-1: automatic convergence detection).)DOC";
constexpr auto _c2py_doc_member_0373262f = R"DOC(Maximum number of warmup cycles when using automatic warmup (safety cap).)DOC";
constexpr auto _c2py_doc_member_96181889 = R"DOC(Seed for random number generator.)DOC";
constexpr auto _c2py_doc_member_c907e317 = R"DOC(Name of random number generator.)DOC";
constexpr auto _c2py_doc_member_01758266 = R"DOC(Maximum runtime in seconds, use -1 to set infinite.)DOC";
constexpr auto _c2py_doc_member_d618aa45 = R"DOC(Verbosity level.)DOC";
constexpr auto _c2py_doc_member_68e6d012 = R"DOC(Add shifting an operator as a move?)DOC";
constexpr auto _c2py_doc_member_3d1a62d7 = R"DOC(Add double insertions as a move?)DOC";
constexpr auto _c2py_doc_member_de609b18 = R"DOC(Calculate the full trace or use an estimate?)DOC";
constexpr auto _c2py_doc_member_19ebb74f = R"DOC(Measure :math:`G(\tau)`? Hermiticity :math:`G_{ij}(\tau) = G_{ji}^*(\tau)` is enforced.)DOC";
constexpr auto _c2py_doc_member_9fe4e2b9 = R"DOC(Measure :math:`G_l` (Legendre)? No hermiticity is enforced.)DOC";
constexpr auto _c2py_doc_member_c478e947 = R"DOC(Measure :math:`O(\tau)` by insertion.)DOC";
constexpr auto _c2py_doc_member_32a7d2a3 = R"DOC(Minimum number of operator insertions in the :math:`O(\tau)` insertion measure.)DOC";
constexpr auto _c2py_doc_member_c35e2239 = R"DOC(Measure :math:`G^{(2)}(\tau,\tau',\tau'')` with three fermionic times.)DOC";
constexpr auto _c2py_doc_member_e514acd9 = R"DOC(Measure :math:`G^{(2)}(i\nu,i\nu',i\nu'')` with three fermionic frequencies.)DOC";
constexpr auto _c2py_doc_member_82de938c = R"DOC(Measure :math:`G^{(2)}(i\nu,i\nu',i\nu'')` with three fermionic frequencies.)DOC";
constexpr auto _c2py_doc_member_94fefae2 = R"DOC(Measure :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-particle channel.)DOC";
constexpr auto _c2py_doc_member_5697b0f5 = R"DOC(Measure :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-particle channel.)DOC";
constexpr auto _c2py_doc_member_8cfeee4a = R"DOC(Measure :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-hole channel.)DOC";
constexpr auto _c2py_doc_member_005cbb6d = R"DOC(Measure :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-hole channel.)DOC";
constexpr auto _c2py_doc_member_6f61d182 = R"DOC(Measure :math:`G^{(2)}(i\omega,l,l')` in the particle-particle channel.)DOC";
constexpr auto _c2py_doc_member_6761c4ea = R"DOC(Measure :math:`G^{(2)}(i\omega,l,l')` in the particle-hole channel.)DOC";
constexpr auto _c2py_doc_member_82215189 = R"DOC(Order of block indices in the definition of :math:`G^{(2)}`.)DOC";
constexpr auto _c2py_doc_member_7d074c93 = R"DOC(List of block index pairs of :math:`G^{(2)}` to measure.)DOC";
constexpr auto _c2py_doc_member_f7409c4c = R"DOC(Number of imaginary-time slices for the :math:`G^{(2)}` measurement.)DOC";
constexpr auto _c2py_doc_member_7a029c0f = R"DOC(Number of bosonic Matsubara frequencies for the :math:`G^{(2)}` measurement.)DOC";
constexpr auto _c2py_doc_member_f0a6b67c = R"DOC(Number of fermionic Matsubara frequencies for the :math:`G^{(2)}` measurement.)DOC";
constexpr auto _c2py_doc_member_e9c77aa2 = R"DOC(Number of Legendre coefficients for the :math:`G^{(2)}(i\omega,l,l')` measurement.)DOC";
constexpr auto _c2py_doc_member_55dd2912 = R"DOC(NFFT buffer size for the :math:`G^{(2)}(i\omega,l,l')` measurement.)DOC";
constexpr auto _c2py_doc_member_27a99f05 = R"DOC(NFFT buffer sizes for different blocks.)DOC";
constexpr auto _c2py_doc_member_e0ebe63f = R"DOC(Measure per-orbital densities via trace-rho-op insertion?)DOC";
constexpr auto _c2py_doc_member_2ba6d7f3 = R"DOC(Measure perturbation order?)DOC";
constexpr auto _c2py_doc_member_bee5d611 = R"DOC(Measure the reduced impurity density matrix?)DOC";
constexpr auto _c2py_doc_member_12ca44c0 = R"DOC(Use the norm of the density matrix in the weight (instead of the trace)?)DOC";
constexpr auto _c2py_doc_member_a2892ac9 = R"DOC(Initial configuration of the run (advanced, use with care).)DOC";
constexpr auto _c2py_doc_member_bb216a99 = R"DOC(Analyse performance of the trace computation with histograms (developers only)?)DOC";
constexpr auto _c2py_doc_member_9c443a14 = R"DOC(Operator insertion/removal probabilities for different blocks.)DOC";
constexpr auto _c2py_doc_member_86410015 =
   R"DOC(List of global moves (with their names). Each move is specified with an index substitution dictionary.)DOC";
constexpr auto _c2py_doc_member_a5b369a5 = R"DOC(Overall probability of the global moves.)DOC";
constexpr auto _c2py_doc_member_3eb4a96e =
   R"DOC(Threshold below which imaginary components of :math:`\Delta` and :math:`h_{loc}` are set to zero.)DOC";
constexpr auto _c2py_doc_member_a4993264 = R"DOC(The maximum size of the determinant matrix before a resize.)DOC";
constexpr auto _c2py_doc_member_0e37c530 = R"DOC(Maximum number of operations before testing the accuracy of :math:`\det(M)` and :math:`M^{-1}`.)DOC";
constexpr auto _c2py_doc_member_e576e2a5 = R"DOC(Threshold for determinant precision warnings.)DOC";
constexpr auto _c2py_doc_member_043bf60b = R"DOC(Threshold for determinant precision error.)DOC";
constexpr auto _c2py_doc_member_559df803 = R"DOC(Bound for the determinant matrix being singular (if :math:`< 0`, checks for subnormal numbers).)DOC";
constexpr auto _c2py_doc_member_569f93b9 = R"DOC(Threshold below which off-diagonal components of :math:`h_{loc}` are set to zero.)DOC";
constexpr auto _c2py_doc_member_5f1bd0b3 =
   R"DOC(Quadratic part of the local Hamiltonian. Must be provided if the :math:`\Delta` interface is used.)DOC";
static PyObject *prop_get_dict_167a5c32(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_a2b1dce8> *)self)->_c);
  c2py::pydict dic;
  dic["h_int"]                         = self_c.h_int;
  dic["n_cycles"]                      = self_c.n_cycles;
  dic["partition_method"]              = self_c.partition_method;
  dic["quantum_numbers"]               = self_c.quantum_numbers;
  dic["loc_n_min"]                     = self_c.loc_n_min;
  dic["loc_n_max"]                     = self_c.loc_n_max;
  dic["length_cycle"]                  = self_c.length_cycle;
  dic["max_length_cycle"]              = self_c.max_length_cycle;
  dic["target_auto_corr_time"]         = self_c.target_auto_corr_time;
  dic["n_warmup_cycles"]               = self_c.n_warmup_cycles;
  dic["max_warmup_cycles"]             = self_c.max_warmup_cycles;
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
  dic["measure_densities"]             = self_c.measure_densities;
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
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_a2b1dce8>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::h_int, _c2py_cls_a2b1dce8>("h_int", _c2py_doc_member_19337922),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::n_cycles, _c2py_cls_a2b1dce8>("n_cycles", _c2py_doc_member_e0206678),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::partition_method, _c2py_cls_a2b1dce8>("partition_method", _c2py_doc_member_2a8cd03e),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::quantum_numbers, _c2py_cls_a2b1dce8>("quantum_numbers", _c2py_doc_member_7ae6201c),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::loc_n_min, _c2py_cls_a2b1dce8>("loc_n_min", _c2py_doc_member_9a5219d0),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::loc_n_max, _c2py_cls_a2b1dce8>("loc_n_max", _c2py_doc_member_903e825a),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::length_cycle, _c2py_cls_a2b1dce8>("length_cycle", _c2py_doc_member_b60dd7e1),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::max_length_cycle, _c2py_cls_a2b1dce8>("max_length_cycle", _c2py_doc_member_d9e01446),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::target_auto_corr_time, _c2py_cls_a2b1dce8>("target_auto_corr_time", _c2py_doc_member_19b5e8ac),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::n_warmup_cycles, _c2py_cls_a2b1dce8>("n_warmup_cycles", _c2py_doc_member_261489c3),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::max_warmup_cycles, _c2py_cls_a2b1dce8>("max_warmup_cycles", _c2py_doc_member_0373262f),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::random_seed, _c2py_cls_a2b1dce8>("random_seed", _c2py_doc_member_96181889),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::random_name, _c2py_cls_a2b1dce8>("random_name", _c2py_doc_member_c907e317),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::max_time, _c2py_cls_a2b1dce8>("max_time", _c2py_doc_member_01758266),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::verbosity, _c2py_cls_a2b1dce8>("verbosity", _c2py_doc_member_d618aa45),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::move_shift, _c2py_cls_a2b1dce8>("move_shift", _c2py_doc_member_68e6d012),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::move_double, _c2py_cls_a2b1dce8>("move_double", _c2py_doc_member_3d1a62d7),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::use_trace_estimator, _c2py_cls_a2b1dce8>("use_trace_estimator", _c2py_doc_member_de609b18),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G_tau, _c2py_cls_a2b1dce8>("measure_G_tau", _c2py_doc_member_19ebb74f),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G_l, _c2py_cls_a2b1dce8>("measure_G_l", _c2py_doc_member_9fe4e2b9),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_O_tau, _c2py_cls_a2b1dce8>("measure_O_tau", _c2py_doc_member_c478e947),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_O_tau_min_ins, _c2py_cls_a2b1dce8>("measure_O_tau_min_ins", _c2py_doc_member_32a7d2a3),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_tau, _c2py_cls_a2b1dce8>("measure_G2_tau", _c2py_doc_member_c35e2239),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_iw, _c2py_cls_a2b1dce8>("measure_G2_iw", _c2py_doc_member_e514acd9),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_iw_nfft, _c2py_cls_a2b1dce8>("measure_G2_iw_nfft", _c2py_doc_member_82de938c),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_iw_pp, _c2py_cls_a2b1dce8>("measure_G2_iw_pp", _c2py_doc_member_94fefae2),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_iw_pp_nfft, _c2py_cls_a2b1dce8>("measure_G2_iw_pp_nfft", _c2py_doc_member_5697b0f5),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_iw_ph, _c2py_cls_a2b1dce8>("measure_G2_iw_ph", _c2py_doc_member_8cfeee4a),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_iw_ph_nfft, _c2py_cls_a2b1dce8>("measure_G2_iw_ph_nfft", _c2py_doc_member_005cbb6d),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_iwll_pp, _c2py_cls_a2b1dce8>("measure_G2_iwll_pp", _c2py_doc_member_6f61d182),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_iwll_ph, _c2py_cls_a2b1dce8>("measure_G2_iwll_ph", _c2py_doc_member_6761c4ea),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_block_order, _c2py_cls_a2b1dce8>("measure_G2_block_order", _c2py_doc_member_82215189),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_blocks, _c2py_cls_a2b1dce8>("measure_G2_blocks", _c2py_doc_member_7d074c93),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_n_tau, _c2py_cls_a2b1dce8>("measure_G2_n_tau", _c2py_doc_member_f7409c4c),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_n_bosonic, _c2py_cls_a2b1dce8>("measure_G2_n_bosonic", _c2py_doc_member_7a029c0f),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_n_fermionic, _c2py_cls_a2b1dce8>("measure_G2_n_fermionic", _c2py_doc_member_f0a6b67c),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_n_l, _c2py_cls_a2b1dce8>("measure_G2_n_l", _c2py_doc_member_e9c77aa2),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_G2_iwll_nfft_buf_size, _c2py_cls_a2b1dce8>("measure_G2_iwll_nfft_buf_size",
                                                                                                       _c2py_doc_member_55dd2912),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::nfft_buf_sizes, _c2py_cls_a2b1dce8>("nfft_buf_sizes", _c2py_doc_member_27a99f05),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_densities, _c2py_cls_a2b1dce8>("measure_densities", _c2py_doc_member_e0ebe63f),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_pert_order, _c2py_cls_a2b1dce8>("measure_pert_order", _c2py_doc_member_2ba6d7f3),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::measure_density_matrix, _c2py_cls_a2b1dce8>("measure_density_matrix", _c2py_doc_member_bee5d611),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::use_norm_as_weight, _c2py_cls_a2b1dce8>("use_norm_as_weight", _c2py_doc_member_12ca44c0),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::initial_configuration, _c2py_cls_a2b1dce8>("initial_configuration", _c2py_doc_member_a2892ac9),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::performance_analysis, _c2py_cls_a2b1dce8>("performance_analysis", _c2py_doc_member_bb216a99),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::proposal_prob, _c2py_cls_a2b1dce8>("proposal_prob", _c2py_doc_member_9c443a14),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::move_global, _c2py_cls_a2b1dce8>("move_global", _c2py_doc_member_86410015),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::move_global_prob, _c2py_cls_a2b1dce8>("move_global_prob", _c2py_doc_member_a5b369a5),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::imag_threshold, _c2py_cls_a2b1dce8>("imag_threshold", _c2py_doc_member_3eb4a96e),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::det_init_size, _c2py_cls_a2b1dce8>("det_init_size", _c2py_doc_member_a4993264),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::det_n_operations_before_check, _c2py_cls_a2b1dce8>("det_n_operations_before_check",
                                                                                                       _c2py_doc_member_0e37c530),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::det_precision_warning, _c2py_cls_a2b1dce8>("det_precision_warning", _c2py_doc_member_e576e2a5),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::det_precision_error, _c2py_cls_a2b1dce8>("det_precision_error", _c2py_doc_member_043bf60b),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::det_singular_threshold, _c2py_cls_a2b1dce8>("det_singular_threshold", _c2py_doc_member_559df803),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::off_diag_threshold, _c2py_cls_a2b1dce8>("off_diag_threshold", _c2py_doc_member_569f93b9),
   c2py::getsetdef_from_member<&_c2py_cls_a2b1dce8::h_loc0, _c2py_cls_a2b1dce8>("h_loc0", _c2py_doc_member_5f1bd0b3),
   {"__dict__", (getter)prop_get_dict_167a5c32, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_a2b1dce8> =
   R"DOC(Parameters passed to the solve method of the solver class.)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_a2b1dce8>;
// --------- class _c2py_cls_7e768e6a -----------
using _c2py_cls_7e768e6a                                            = triqs_cthyb::solver_core;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_7e768e6a>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_7e768e6a> = "triqs_cthyb.solver_core.SolverCore";
static const auto _c2py_init_70a39a57 =
   c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_7e768e6a, const triqs_cthyb::constr_parameters_t &>("p")};
template <> constexpr initproc c2py::tp_init<_c2py_cls_7e768e6a> = c2py::pyfkw_constructor<_c2py_init_70a39a57>;
template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_7e768e6a> =
   _c2py_init_70a39a57.doc(R"DOC(
Construct a CTHYB solver.

Parameters
----------
p : {par_0}
   Parameters used for constructing the solver.
)DOC",
                           {{c2py::python_typename<const triqs_cthyb::constr_parameters_t &>()}});
// solve
static auto const _c2py_fun_cdee5b11 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_7e768e6a &self, const triqs_cthyb::solve_parameters_t &p) -> decltype(auto) { return self.solve(p); }, "self", "p")};

static const auto _c2py_doc_cdee5b11 = _c2py_fun_cdee5b11.doc(R"DOC(
Solve the impurity problem.

Parameters
----------
p : {par_0}
   Parameters controlling the Monte Carlo simulation and measurements.
)DOC",
                                                              {{c2py::python_typename<const triqs_cthyb::solve_parameters_t &>()}});

// ----- Method table ----
// clang-format off
template <> PyMethodDef c2py::tp_methods<_c2py_cls_7e768e6a>[] = {
   PMDF("solve", cdee5b11),
   {"__write_hdf5__", c2py::tpxx_write_h5<_c2py_cls_7e768e6a>, METH_VARARGS, "  "},
   {"__getstate__", c2py::getstate_h5<_c2py_cls_7e768e6a>, METH_NOARGS, ""},
   {"__setstate__", c2py::setstate_h5<_c2py_cls_7e768e6a>, METH_O, ""},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};
// clang-format on

constexpr auto _c2py_doc_member_28c40d2e = R"DOC(Parameters used for constructing the solver.)DOC";
constexpr auto _c2py_doc_member_6ca31b4a = R"DOC(Parameters passed to the solve method.)DOC";
constexpr auto _c2py_doc_member_a630495a = R"DOC(Single-particle Green's function :math:`G(\tau)` in imaginary time.)DOC";
constexpr auto _c2py_doc_member_789d035e = R"DOC(Intermediate Green's function used to accumulate :math:`G(\tau)` (real or complex).)DOC";
constexpr auto _c2py_doc_member_c4b4ab9c = R"DOC(Violation of the property :math:`G_{ij}(\tau) = G_{ji}^*(\tau)` after the measurement.)DOC";
constexpr auto _c2py_doc_member_bb9bd800 = R"DOC(Single-particle Green's function :math:`G_l` in the Legendre representation.)DOC";
constexpr auto _c2py_doc_member_2fc3b9f2 = R"DOC(General operator Green's function :math:`O(\tau)` in imaginary time.)DOC";
constexpr auto _c2py_doc_member_7651e12a = R"DOC(Two-particle Green's function :math:`G^{(2)}(\tau_1,\tau_2,\tau_3)` with three fermionic times.)DOC";
constexpr auto _c2py_doc_member_411ac008 =
   R"DOC(Two-particle Green's function :math:`G^{(2)}(i\nu,i\nu',i\nu'')` with three fermionic frequencies.)DOC";
constexpr auto _c2py_doc_member_2e209577 =
   R"DOC(Two-particle Green's function :math:`G^{(2)}(i\nu,i\nu',i\nu'')` with three fermionic frequencies.)DOC";
constexpr auto _c2py_doc_member_d66754fd =
   R"DOC(Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-particle channel.)DOC";
constexpr auto _c2py_doc_member_d922ac20 =
   R"DOC(Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-particle channel.)DOC";
constexpr auto _c2py_doc_member_de676195 = R"DOC(Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-hole channel.)DOC";
constexpr auto _c2py_doc_member_f5f6f4a8 = R"DOC(Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the particle-hole channel.)DOC";
constexpr auto _c2py_doc_member_f1ffaa95 = R"DOC(Two-particle Green's function :math:`G^{(2)}(i\omega,l,l')` in the particle-particle channel.)DOC";
constexpr auto _c2py_doc_member_e9ff9dfd = R"DOC(Two-particle Green's function :math:`G^{(2)}(i\omega,l,l')` in the particle-hole channel.)DOC";
constexpr auto _c2py_doc_member_5bfd95dd = R"DOC(Histogram of the total perturbation order.)DOC";
constexpr auto _c2py_doc_member_b85d7f38 = R"DOC(Histograms of the perturbation order for each block.)DOC";
static constexpr auto prop_doc_3bddb197  = R"DOC(:math:`G_0^{-1}(i\omega_n = \infty)` in Matsubara frequencies.)DOC";
static constexpr auto prop_doc_8edb23e9  = R"DOC(Hybridization function :math:`\Delta(\tau)` in imaginary time.)DOC";
static constexpr auto prop_doc_2748acea  = R"DOC(Non-interacting Green's function :math:`G_0(i\omega)` in Matsubara frequencies.)DOC";
static constexpr auto prop_doc_2036d594  = R"DOC(Auto-correlation time in units of MC cycles.)DOC";
static constexpr auto prop_doc_4545a854 =
   R"DOC(Whether the auto-correlation time estimate has saturated (false: it is only a lower bound, run longer).)DOC";
static constexpr auto prop_doc_d6d80776 = R"DOC(Average perturbation order.)DOC";
static constexpr auto prop_doc_41d238d3 = R"DOC(Error bar for average perturbation order.)DOC";
static constexpr auto prop_doc_45b52fcb = R"DOC(Monte Carlo average sign.)DOC";
static constexpr auto prop_doc_6d131e0a = R"DOC(Error bar for average sign.)DOC";
static constexpr auto prop_doc_12d5fb8a = R"DOC(Per-orbital densities, organized by blocks)DOC";
static constexpr auto prop_doc_79de08a4 = R"DOC(Error bars for per-orbital densities, organized by blocks)DOC";
static constexpr auto prop_doc_cfd81bc2 = R"DOC(Accumulated density matrix.)DOC";
static constexpr auto prop_doc_0650144c = R"DOC(Error bars for density matrix)DOC";
static constexpr auto prop_doc_dd280cdb = R"DOC(The local Hamiltonian :math:`H_{loc}` used in the last solve.)DOC";
static constexpr auto prop_doc_110c55f1 = R"DOC(The noninteracting part of the local Hamiltonian.)DOC";
static constexpr auto prop_doc_e9d8c699 = R"DOC(Diagonalization of :math:`H_{loc}`.)DOC";
static constexpr auto prop_doc_5384c449 = R"DOC(Is the solver compiled with support for complex hybridization?)DOC";
static constexpr auto prop_doc_7254e6c7 = R"DOC(Final configuration of the last solve call.)DOC";
static constexpr auto prop_doc_af94c241 = R"DOC(Parameters used for constructing the solver.)DOC";
static constexpr auto prop_doc_92a8d1b7 = R"DOC(Parameters used in the last solve.)DOC";
static constexpr auto prop_doc_a78e6b75 = R"DOC(Effective length_cycle used in the accumulation (may differ from length_cycle in auto mode).)DOC";
static constexpr auto prop_doc_31d56a88 = R"DOC(Is the solver compiled with support for a complex local Hamiltonian?)DOC";
static constexpr auto prop_doc_5aaceaa3 = R"DOC(Histograms related to the performance analysis.)DOC";
static constexpr auto prop_doc_33aa8ad4 = R"DOC(Status of the ``solve()`` on exit.)DOC";
static constexpr auto prop_doc_3daf94b1 = R"DOC(Number of warmup cycles actually performed (may differ from n_warmup_cycles in auto mode).)DOC";

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_7e768e6a>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::constr_parameters, _c2py_cls_7e768e6a>("constr_parameters", _c2py_doc_member_28c40d2e),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::solve_parameters, _c2py_cls_7e768e6a>("solve_parameters", _c2py_doc_member_6ca31b4a),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::G_tau, _c2py_cls_7e768e6a>("G_tau", _c2py_doc_member_a630495a),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::G_tau_accum, _c2py_cls_7e768e6a>("G_tau_accum", _c2py_doc_member_789d035e),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::asymmetry_G_tau, _c2py_cls_7e768e6a>("asymmetry_G_tau", _c2py_doc_member_c4b4ab9c),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::G_l, _c2py_cls_7e768e6a>("G_l", _c2py_doc_member_bb9bd800),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::O_tau, _c2py_cls_7e768e6a>("O_tau", _c2py_doc_member_2fc3b9f2),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::G2_tau, _c2py_cls_7e768e6a>("G2_tau", _c2py_doc_member_7651e12a),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::G2_iw, _c2py_cls_7e768e6a>("G2_iw", _c2py_doc_member_411ac008),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::G2_iw_nfft, _c2py_cls_7e768e6a>("G2_iw_nfft", _c2py_doc_member_2e209577),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::G2_iw_pp, _c2py_cls_7e768e6a>("G2_iw_pp", _c2py_doc_member_d66754fd),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::G2_iw_pp_nfft, _c2py_cls_7e768e6a>("G2_iw_pp_nfft", _c2py_doc_member_d922ac20),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::G2_iw_ph, _c2py_cls_7e768e6a>("G2_iw_ph", _c2py_doc_member_de676195),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::G2_iw_ph_nfft, _c2py_cls_7e768e6a>("G2_iw_ph_nfft", _c2py_doc_member_f5f6f4a8),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::G2_iwll_pp, _c2py_cls_7e768e6a>("G2_iwll_pp", _c2py_doc_member_f1ffaa95),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::G2_iwll_ph, _c2py_cls_7e768e6a>("G2_iwll_ph", _c2py_doc_member_e9ff9dfd),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::perturbation_order_total, _c2py_cls_7e768e6a>("perturbation_order_total",
                                                                                                  _c2py_doc_member_5bfd95dd),
   c2py::getsetdef_from_member<&_c2py_cls_7e768e6a::perturbation_order, _c2py_cls_7e768e6a>("perturbation_order", _c2py_doc_member_b85d7f38),
   {"Delta_infty", c2py::getter_from_method<c2py::castm<>(&triqs_cthyb::solver_core::Delta_infty)>, nullptr, prop_doc_3bddb197, nullptr},
   {"Delta_tau", c2py::getter_from_method<c2py::castm<>(&triqs_cthyb::solver_core::Delta_tau)>, nullptr, prop_doc_8edb23e9, nullptr},
   {"G0_iw", c2py::getter_from_method<c2py::castm<>(&triqs_cthyb::solver_core::G0_iw)>, nullptr, prop_doc_2748acea, nullptr},
   {"auto_corr_time", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::auto_corr_time)>, nullptr, prop_doc_2036d594, nullptr},
   {"auto_corr_time_converged", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::auto_corr_time_converged)>, nullptr,
    prop_doc_4545a854, nullptr},
   {"average_order", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::average_order)>, nullptr, prop_doc_d6d80776, nullptr},
   {"average_order_error", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::average_order_error)>, nullptr, prop_doc_41d238d3,
    nullptr},
   {"average_sign", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::average_sign)>, nullptr, prop_doc_45b52fcb, nullptr},
   {"average_sign_error", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::average_sign_error)>, nullptr, prop_doc_6d131e0a,
    nullptr},
   {"densities", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::densities)>, nullptr, prop_doc_12d5fb8a, nullptr},
   {"densities_errors", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::densities_errors)>, nullptr, prop_doc_79de08a4, nullptr},
   {"density_matrix", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::density_matrix)>, nullptr, prop_doc_cfd81bc2, nullptr},
   {"density_matrix_errors", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::density_matrix_errors)>, nullptr, prop_doc_0650144c,
    nullptr},
   {"h_loc", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::h_loc)>, nullptr, prop_doc_dd280cdb, nullptr},
   {"h_loc0", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::h_loc0)>, nullptr, prop_doc_110c55f1, nullptr},
   {"h_loc_diagonalization", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::h_loc_diagonalization)>, nullptr, prop_doc_e9d8c699,
    nullptr},
   {"hybridisation_is_complex", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::hybridisation_is_complex)>, nullptr,
    prop_doc_5384c449, nullptr},
   {"last_configuration", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::last_configuration)>, nullptr, prop_doc_7254e6c7,
    nullptr},
   {"last_constr_parameters", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::last_constr_parameters)>, nullptr, prop_doc_af94c241,
    nullptr},
   {"last_solve_parameters", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::last_solve_parameters)>, nullptr, prop_doc_92a8d1b7,
    nullptr},
   {"length_cycle_used", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::length_cycle_used)>, nullptr, prop_doc_a78e6b75, nullptr},
   {"local_hamiltonian_is_complex", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::local_hamiltonian_is_complex)>, nullptr,
    prop_doc_31d56a88, nullptr},
   {"performance_analysis", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::get_performance_analysis)>, nullptr, prop_doc_5aaceaa3,
    nullptr},
   {"solve_status", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::solve_status)>, nullptr, prop_doc_33aa8ad4, nullptr},
   {"warmup_cycles_done", c2py::getter_from_method<c2py::castmc<>(&triqs_cthyb::solver_core::warmup_cycles_done)>, nullptr, prop_doc_3daf94b1,
    nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_7e768e6a> = R"DOC(Continuous-time hybridization-expansion quantum Monte Carlo solver.)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_7e768e6a>;

// ==================== module functions ====================

//--------------------- module function table  -----------------------------

// clang-format off
static PyMethodDef module_methods[] = {
   {nullptr, nullptr, 0, nullptr}  // Sentinel
};
// clang-format on

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

  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_e57534ae>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_a2b1dce8>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_7e768e6a>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  if (not c2py::register_internal_types()) return NULL;
#define _add_type(T, N)                                                                                                                              \
  if (not c2py::add_type_object_to_main<T>(N, m)) return NULL
  _add_type(_c2py_cls_e57534ae, "ConstrParametersT");
  _add_type(_c2py_cls_a2b1dce8, "SolveParametersT");
  _add_type(_c2py_cls_7e768e6a, "SolverCore");
#undef _add_type

  c2py::pyref module = c2py::pyref::module("h5.formats");
  if (not module) return nullptr;
  c2py::pyref register_class = module.attr("register_class");

  register_h5_type<_c2py_cls_7e768e6a>(register_class);

  return m;
}
#endif
// CLAIR_WRAP_GEN
