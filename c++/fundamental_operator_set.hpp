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
#ifndef TRIQS_CTQMC_KRYLOV_FUNDAMENTAL_OPERATOR_SET
#define TRIQS_CTQMC_KRYLOV_FUNDAMENTAL_OPERATOR_SET
//#include <utility>
#include <vector>
#include <map>
#include "./operator.hpp"

namespace cthyb_krylov {

// This class defines the list of operators that are used to describe
// e.g. Fock states etc.

class fundamental_operator_set {
 public:
 using indices_t = triqs::utility::many_body_operator<double>::indices_t;

 private:
 // the table index <-> n
 using map_t = std::map<indices_t, int>;
 map_t map_index_n;

 public:
 fundamental_operator_set() {}

 // constructor on a vector which gives the number of alpha's for every a
 // this only makes sense if the indices are two integers
 fundamental_operator_set(std::vector<int> const& v) {
  for (int i = 0; i < v.size(); ++i)
   for (int j = 0; j < v[i]; ++j) add_operator(i, j);
 }

 // REMOVE THIS : jsut in construction, to avoid reorganizing all the time
 // may lead to bug if start to use it, then add new ops, then reuse...
 // add an operator in the set
 template <typename... IndexType> void add_operator(IndexType const&... ind) {
  map_index_n.insert({{ind...}, n_operators()});
  // reorder the indices which are always given in the order of the indices tuple
  map_t m;
  int i = 0;
  for (auto const& p : map_index_n) m.insert({p.first, i++});
  std::swap(m, map_index_n);
 }

 // return the dimension of the space spanned by the operators
 int dimension() const { return 1ull << n_operators(); } // 2^ n_ops

 // return the number of operators
 int n_operators() const { return map_index_n.size(); }

 // flatten (a,alpha) --> n
 int operator[](indices_t const& t) const { return map_index_n.at(t); }

 // iterator on the tuples
 // for (auto & x : fops) { x.linear_index is linear_index, while x.index is the C multi-index.
 struct _cdress {
  int linear_index;
  indices_t const& index;
  _cdress(typename map_t::const_iterator _it) : index(_it->first), linear_index(_it->second) {}
 };
 using const_iterator = triqs::utility::dressed_iterator<typename map_t::const_iterator, _cdress>;

 const_iterator begin() const noexcept { return map_index_n.begin(); }
 const_iterator end() const noexcept { return map_index_n.end(); }
 const_iterator cbegin() const noexcept { return map_index_n.cbegin(); }
 const_iterator cend() const noexcept { return map_index_n.cend(); }
};
}
#endif
