/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by I. Krivenko, M. Ferrero, O. Parcollet
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

#ifndef TRIQS_CTQMC_KRYLOV_TUPLE_HASHER
#define TRIQS_CTQMC_KRYLOV_TUPLE_HASHER

#include <tuple>
#include <boost/functional/hash.hpp>
#include <triqs/utility/tuple_tools.hpp>

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

// because there is no std hash for things like tuple<int,int>
template<typename TupleType>
struct tuple_hasher: std::unary_function<TupleType, size_t> {

  struct combine_wrap {
    combine_wrap(size_t s_): s(s_) {}
    template<typename T>
    void operator() (T x) { boost::hash_combine(s, x); }
    size_t s;
  };

  size_t operator()(TupleType const & t) const {
    size_t seed = 0;
    triqs::tuple::for_each(t, combine_wrap(seed));
    return seed;
  }

};

}}}}
#endif
