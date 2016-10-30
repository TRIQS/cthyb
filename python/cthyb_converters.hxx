// DO NOT EDIT
// Generated automatically using libclang using the command :
// c++2py.py ../c++/solver_core.hpp -I../../cthyb.build/c++ -I../c++ -p -mpytriqs.applications.impurity_solvers.cthyb -o cthyb --moduledoc "The cthyb solver"


// --- C++ Python converter for solve_parameters_t
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/string.hpp>
#include <algorithm>

namespace triqs { namespace py_tools {

template <> struct py_converter<solve_parameters_t> {
 static PyObject *c2py(solve_parameters_t const & x) {
  PyObject * d = PyDict_New();
  PyDict_SetItemString( d, "h_int"                 , convert_to_python(x.h_int));
  PyDict_SetItemString( d, "n_cycles"              , convert_to_python(x.n_cycles));
  PyDict_SetItemString( d, "partition_method"      , convert_to_python(x.partition_method));
  PyDict_SetItemString( d, "quantum_numbers"       , convert_to_python(x.quantum_numbers));
  PyDict_SetItemString( d, "length_cycle"          , convert_to_python(x.length_cycle));
  PyDict_SetItemString( d, "n_warmup_cycles"       , convert_to_python(x.n_warmup_cycles));
  PyDict_SetItemString( d, "random_seed"           , convert_to_python(x.random_seed));
  PyDict_SetItemString( d, "random_name"           , convert_to_python(x.random_name));
  PyDict_SetItemString( d, "max_time"              , convert_to_python(x.max_time));
  PyDict_SetItemString( d, "verbosity"             , convert_to_python(x.verbosity));
  PyDict_SetItemString( d, "move_shift"            , convert_to_python(x.move_shift));
  PyDict_SetItemString( d, "move_double"           , convert_to_python(x.move_double));
  PyDict_SetItemString( d, "use_trace_estimator"   , convert_to_python(x.use_trace_estimator));
  PyDict_SetItemString( d, "measure_g_tau"         , convert_to_python(x.measure_g_tau));
  PyDict_SetItemString( d, "measure_g_l"           , convert_to_python(x.measure_g_l));
  PyDict_SetItemString( d, "measure_pert_order"    , convert_to_python(x.measure_pert_order));
  PyDict_SetItemString( d, "measure_density_matrix", convert_to_python(x.measure_density_matrix));
  PyDict_SetItemString( d, "use_norm_as_weight"    , convert_to_python(x.use_norm_as_weight));
  PyDict_SetItemString( d, "performance_analysis"  , convert_to_python(x.performance_analysis));
  PyDict_SetItemString( d, "proposal_prob"         , convert_to_python(x.proposal_prob));
  PyDict_SetItemString( d, "move_global"           , convert_to_python(x.move_global));
  PyDict_SetItemString( d, "move_global_prob"      , convert_to_python(x.move_global_prob));
  PyDict_SetItemString( d, "imag_threshold"        , convert_to_python(x.imag_threshold));
  return d;
 }

 template <typename T, typename U> static void _get_optional(PyObject *dic, const char *name, T &r, U const &init_default) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = init_default;
 }

 template <typename T> static void _get_optional(PyObject *dic, const char *name, T &r) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = T{};
 }

 static solve_parameters_t py2c(PyObject *dic) {
  solve_parameters_t res;
  res.h_int = convert_from_python<many_body_op_t>(PyDict_GetItemString(dic, "h_int"));
  res.n_cycles = convert_from_python<int>(PyDict_GetItemString(dic, "n_cycles"));
  _get_optional(dic, "partition_method"      , res.partition_method         ,"autopartition");
  _get_optional(dic, "quantum_numbers"       , res.quantum_numbers          ,std::vector<many_body_op_t>{});
  _get_optional(dic, "length_cycle"          , res.length_cycle             ,50);
  _get_optional(dic, "n_warmup_cycles"       , res.n_warmup_cycles          ,5000);
  _get_optional(dic, "random_seed"           , res.random_seed              ,34788+928374*triqs::mpi::communicator().rank());
  _get_optional(dic, "random_name"           , res.random_name              ,"");
  _get_optional(dic, "max_time"              , res.max_time                 ,-1);
  _get_optional(dic, "verbosity"             , res.verbosity                ,((triqs::mpi::communicator().rank()==0)?3:0));
  _get_optional(dic, "move_shift"            , res.move_shift               ,true);
  _get_optional(dic, "move_double"           , res.move_double              ,false);
  _get_optional(dic, "use_trace_estimator"   , res.use_trace_estimator      ,false);
  _get_optional(dic, "measure_g_tau"         , res.measure_g_tau            ,true);
  _get_optional(dic, "measure_g_l"           , res.measure_g_l              ,false);
  _get_optional(dic, "measure_pert_order"    , res.measure_pert_order       ,false);
  _get_optional(dic, "measure_density_matrix", res.measure_density_matrix   ,false);
  _get_optional(dic, "use_norm_as_weight"    , res.use_norm_as_weight       ,false);
  _get_optional(dic, "performance_analysis"  , res.performance_analysis     ,false);
  _get_optional(dic, "proposal_prob"         , res.proposal_prob            ,(std::map<std::string,double>{}));
  _get_optional(dic, "move_global"           , res.move_global              ,(std::map<std::string,indices_map_t>{}));
  _get_optional(dic, "move_global_prob"      , res.move_global_prob         ,0.05);
  _get_optional(dic, "imag_threshold"        , res.imag_threshold           ,1.e-15);
  return res;
 }

 template <typename T>
 static void _check(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!convertible_from_python<T>(PyDict_GetItemString(dic, name), false))
   fs << "\n" << ++err << " The parameter " << name << " does not have the right type : expecting " << tname
      << " in C++, but got '" << PyDict_GetItemString(dic, name)->ob_type->tp_name << "' in Python.";
 }

 template <typename T>
 static void _check_mandatory(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!PyDict_Contains(dic, pyref::string(name)))
   fs << "\n" << ++err << " Mandatory parameter " << name << " is missing.";
  else _check<T>(dic,fs,err,name,tname);
 }

 template <typename T>
 static void _check_optional(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (PyDict_Contains(dic, pyref::string(name))) _check<T>(dic, fs, err, name, tname);
 }

 static bool is_convertible(PyObject *dic, bool raise_exception) {
  if (dic == nullptr or !PyDict_Check(dic)) {
   if (raise_exception) { PyErr_SetString(PyExc_TypeError, "The function must be called with named arguments");}
   return false;
  }
  std::stringstream fs, fs2; int err=0;

#ifndef TRIQS_ALLOW_UNUSED_PARAMETERS
  std::vector<std::string> ks, all_keys = {"h_int","n_cycles","partition_method","quantum_numbers","length_cycle","n_warmup_cycles","random_seed","random_name","max_time","verbosity","move_shift","move_double","use_trace_estimator","measure_g_tau","measure_g_l","measure_pert_order","measure_density_matrix","use_norm_as_weight","performance_analysis","proposal_prob","move_global","move_global_prob","imag_threshold"};
  pyref keys = PyDict_Keys(dic);
  if (!convertible_from_python<std::vector<std::string>>(keys, true)) {
   fs << "\nThe dict keys are not strings";
   goto _error;
  }
  ks = convert_from_python<std::vector<std::string>>(keys);
  for (auto & k : ks)
   if (std::find(all_keys.begin(), all_keys.end(), k) == all_keys.end())
    fs << "\n"<< ++err << " The parameter '" << k << "' is not recognized.";
#endif

  _check_mandatory<many_body_op_t                      >(dic, fs, err, "h_int"                 , "many_body_op_t");
  _check_mandatory<int                                 >(dic, fs, err, "n_cycles"              , "int");
  _check_optional <std::string                         >(dic, fs, err, "partition_method"      , "std::string");
  _check_optional <std::vector<many_body_op_t>         >(dic, fs, err, "quantum_numbers"       , "std::vector<many_body_op_t>");
  _check_optional <int                                 >(dic, fs, err, "length_cycle"          , "int");
  _check_optional <int                                 >(dic, fs, err, "n_warmup_cycles"       , "int");
  _check_optional <int                                 >(dic, fs, err, "random_seed"           , "int");
  _check_optional <std::string                         >(dic, fs, err, "random_name"           , "std::string");
  _check_optional <int                                 >(dic, fs, err, "max_time"              , "int");
  _check_optional <int                                 >(dic, fs, err, "verbosity"             , "int");
  _check_optional <bool                                >(dic, fs, err, "move_shift"            , "bool");
  _check_optional <bool                                >(dic, fs, err, "move_double"           , "bool");
  _check_optional <bool                                >(dic, fs, err, "use_trace_estimator"   , "bool");
  _check_optional <bool                                >(dic, fs, err, "measure_g_tau"         , "bool");
  _check_optional <bool                                >(dic, fs, err, "measure_g_l"           , "bool");
  _check_optional <bool                                >(dic, fs, err, "measure_pert_order"    , "bool");
  _check_optional <bool                                >(dic, fs, err, "measure_density_matrix", "bool");
  _check_optional <bool                                >(dic, fs, err, "use_norm_as_weight"    , "bool");
  _check_optional <bool                                >(dic, fs, err, "performance_analysis"  , "bool");
  _check_optional <std::map<std::string, double>       >(dic, fs, err, "proposal_prob"         , "std::map<std::string, double>");
  _check_optional <std::map<std::string, indices_map_t>>(dic, fs, err, "move_global"           , "std::map<std::string, indices_map_t>");
  _check_optional <double                              >(dic, fs, err, "move_global_prob"      , "double");
  _check_optional <double                              >(dic, fs, err, "imag_threshold"        , "double");
  if (err) goto _error;
  return true;

 _error:
   fs2 << "\n---- There " << (err > 1 ? "are " : "is ") << err<< " error"<<(err >1 ?"s" : "")<< " in Python -> C++ transcription for the class solve_parameters_t\n" <<fs.str();
   if (raise_exception) PyErr_SetString(PyExc_TypeError, fs2.str().c_str());
  return false;
 }
};

}}