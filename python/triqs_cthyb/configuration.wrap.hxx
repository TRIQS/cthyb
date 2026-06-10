#include <c2py/c2py.hpp>

#ifndef C2PY_HXX_DECLARATION_configuration_GUARDS
#define C2PY_HXX_DECLARATION_configuration_GUARDS
template <> constexpr bool c2py::is_wrapped<triqs_cthyb::op_desc>           = true;
template <> inline constexpr auto c2py::tp_name<triqs_cthyb::op_desc>       = "triqs_cthyb.configuration.OpDesc";
template <> constexpr bool c2py::is_wrapped<triqs_cthyb::configuration>     = true;
template <> inline constexpr auto c2py::tp_name<triqs_cthyb::configuration> = "triqs_cthyb.configuration.Configuration";
#endif