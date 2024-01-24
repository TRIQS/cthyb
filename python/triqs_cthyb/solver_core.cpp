#include <c2py/c2py.hpp>
#include "triqs_cthyb/solver_core.hpp"
#include <triqs/stat.hpp>
#include <triqs/operators.hpp>
#include <triqs/atom_diag.hpp>

namespace c2py_module {

  auto documentation = "The TRIQS cthyb solver";
  auto package_name  = "triqs_cthyb";

  auto get_set_as_properties = true;

  auto match_names =  "triqs_cthyb::(solver_core|solve_parameters_t|constr_parameters_t|block_order)";


} // namespace c2py_module
