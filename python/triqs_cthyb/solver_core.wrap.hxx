#include <c2py/c2py.hpp>

#ifndef C2PY_HXX_DECLARATION_solver_core_GUARDS
#define C2PY_HXX_DECLARATION_solver_core_GUARDS
template <> constexpr bool c2py::is_wrapped<triqs_cthyb::constr_parameters_t>     = true;
template <> inline constexpr auto c2py::tp_name<triqs_cthyb::constr_parameters_t> = "triqs_cthyb.solver_core.ConstrParametersT";
template <> constexpr bool c2py::is_wrapped<triqs_cthyb::solve_parameters_t>      = true;
template <> inline constexpr auto c2py::tp_name<triqs_cthyb::solve_parameters_t>  = "triqs_cthyb.solver_core.SolveParametersT";
template <> constexpr bool c2py::is_wrapped<triqs_cthyb::solver_core>             = true;
template <> inline constexpr auto c2py::tp_name<triqs_cthyb::solver_core>         = "triqs_cthyb.solver_core.SolverCore";
template <> constexpr bool c2py::is_wrapped<triqs_cthyb::block_order>             = true;
template <>
const std::map<triqs_cthyb::block_order, str_t> c2py::enum_to_string<triqs_cthyb::block_order> = {{triqs_cthyb::block_order::AABB, "AABB"},
                                                                                                  {triqs_cthyb::block_order::ABBA, "ABBA"}};
#endif