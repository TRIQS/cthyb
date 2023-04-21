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

  // FIXME : error in fun_accept clause...
  auto fun_accept_only  = "solver_core";
  auto enum_accept_only = "triqs_cthyb::block_order";

  template <> struct wrap_info<triqs_cthyb::solver_core> {
    using bases_to_merge = std::tuple<triqs_cthyb::container_set_t>;
  };

  template <> struct wrap_info<triqs_cthyb::constr_parameters_t> {
    static constexpr auto synthetize_init_from_pydict = true;
  };

  template <> struct wrap_info<triqs_cthyb::solve_parameters_t> {
    static constexpr auto synthetize_init_from_pydict = true;
  };

} // namespace c2py_module
