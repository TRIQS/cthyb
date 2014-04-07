/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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
#include "gf_block_structure.hpp"

namespace cthyb_matrix {

gf_block_structure_t::gf_block_structure_t(fundamental_operator_set const& fops, std::vector<block_desc_t> const& block_structure) {

 using indices_t = fundamental_operator_set::indices_t;
 std::map<indices_t, std::pair<int, int>> indices_to_ints;
 for (int bl = 0; bl < block_structure.size(); ++bl) {
  auto const& indices = block_structure[bl].indices;
  for (int i = 0; i < indices.size(); ++i) {
   indices_to_ints[indices[i]] = std::make_pair(bl, i);
  }
 }

 // create the map int_pair_to_n : (int,int) --> int identifying operators
 for (auto x : fops) int_pair_to_n[indices_to_ints.at(x.index)] = x.linear_index;
}
}
