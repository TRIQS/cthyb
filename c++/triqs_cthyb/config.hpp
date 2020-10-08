/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once
#include <triqs/arrays.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <string>
#include <vector>
#include <map>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/state.hpp>
#include <triqs/hilbert_space/imperative_operator.hpp>
//#include "./array_suppl.hpp"

#include <h5/h5.hpp>

namespace triqs_cthyb {

using namespace triqs::arrays;
using namespace triqs::gfs;
using namespace triqs::mesh;
using namespace triqs::hilbert_space;
namespace operators = triqs::operators;

using std::isfinite;
inline bool isfinite(dcomplex const &x) { return std::isfinite(real(x)) && std::isfinite(imag(x)); }

inline double real(double x) { return x; }
inline double imag(double) { return 0; }

// FIXME when moved into the lib, we have to template of h_scalar_t = double and complex
// consider moving this into the class, under the template

#ifdef HYBRIDISATION_IS_COMPLEX
using det_scalar_t = dcomplex;
using delta_target_t = matrix_valued;
#else
using det_scalar_t = double;
using delta_target_t = matrix_real_valued;
#endif

#ifdef LOCAL_HAMILTONIAN_IS_COMPLEX
using h_scalar_t = dcomplex; // type of scalar for H_loc: double or complex.
static constexpr bool is_h_scalar_complex = true;
#else
using h_scalar_t = double; // type of scalar for H_loc: double or complex.
static constexpr bool is_h_scalar_complex = false;
#endif

using mc_weight_t = decltype(h_scalar_t{} * det_scalar_t{}); // complex iif either is complex
using many_body_op_t = triqs::operators::many_body_operator_generic<h_scalar_t>; // Operator with real or complex value
using matrix_t = matrix<h_scalar_t>;
using G_target_t = std::conditional_t<triqs::is_complex<mc_weight_t>::value, matrix_valued, matrix_real_valued>;

}
