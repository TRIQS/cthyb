#include <c2py/c2py.hpp>

#include <triqs/operators.hpp>
#include <triqs/atom_diag.hpp>
#include <triqs/stat.hpp>
#include <triqs/mesh.hpp>
#include <triqs/gfs.hpp>
#include <triqs/c2py_converters/gf.hpp>
#include <triqs/c2py_converters/mesh.hpp>
#include <triqs/c2py_converters/fundamental_operator_set.hpp>
#include <triqs/c2py_converters/real_or_complex.hpp>
#include <triqs/c2py_converters/operators_real_complex.hpp>
#include <triqs/c2py_converters/arrays.hpp>
#include <triqs_cthyb/solver_core.hpp>

#include "./configuration.wrap.hxx"

#include "solver_core.wrap.cxx"
