#include "triqs_cthyb/solver_core.hpp"

#include <c2py/c2py.hpp>

// ------------  ALL GOES INTO TRIQS _-----------------
// FIXME auto include in TRIQS under macro CLAIR_C2PY_MODULE
// define automatically in compiler pass ?
#include <nda_py/cpp2py_converters.hpp>
#include <triqs/cpp2py_converters.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>

// FIXME : should be in TRIQS itself ... Same macro
template <> static constexpr bool c2py::is_wrapped<triqs::operators::many_body_operator>                          = true;
template <> static constexpr bool c2py::is_wrapped<triqs::stat::histogram>                                        = true;
template <> static constexpr bool c2py::is_wrapped<triqs::atom_diag::atom_diag<triqs_cthyb::is_h_scalar_complex>> = true;
// ------------  END ALL GOES INTO TRIQS _-----------------


namespace c2py_module {

  auto documentation = "The TRIQS cthyb solver";
  auto ns            = "triqs_cthyb";
  auto package_name  = "triqs_cthyb";

  auto get_set_as_properties = true;

  auto regex_filter =  "triqs_cthyb::(solver_core|solve_parameters_t|constr_parameters_t|block_order)";

  // Alternatively, one can be more specific fun/cls/enum... Same result
  //auto regex_exclude_fun  = ".*"; 
  //auto regex_filter_cls  = "triqs_cthyb::(solver_core|solve_parameters_t|constr_parameters_t)";
  //auto regex_filter_enum = "triqs_cthyb::block_order";


} // namespace c2py_module
