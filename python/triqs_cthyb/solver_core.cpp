#include <c2py/c2py.hpp>
#include "triqs_cthyb/solver_core.hpp"
#include <triqs/stat.hpp>
#include <triqs/operators.hpp>
#include <triqs/atom_diag.hpp>

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
