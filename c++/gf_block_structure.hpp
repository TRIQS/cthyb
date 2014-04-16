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
#pragma once
#include <triqs/arrays.hpp>
#include <triqs/gfs.hpp>
#include <string>
#include <vector>
#include <map>

#include "triqs/draft/hilbert_space_tools/fundamental_operator_set.hpp"

using namespace triqs::arrays;

namespace cthyb {

/// ???
struct block_desc_t {
 std::string name;
 std::vector<fundamental_operator_set::indices_t> indices;
 block_desc_t(std::string name, std::vector<fundamental_operator_set::indices_t> ind) : name(name), indices(ind){}
 block_desc_t(std::string name="") : name(name), indices(){}
 // for python interface only
 void indices_push_back(std::string a, std::string b) {
  indices.push_back({a, b});
 }
};

struct gf_block_structure_t { 

 // Constructor
 gf_block_structure_t(fundamental_operator_set const& fops, std::vector<block_desc_t> const& block_structure);

 // get fundamental operators
 long get_fundamental_operator_linear_index(int block_index, int inner_index) const {
  return int_pair_to_n.at({block_index, inner_index});
 }

 private:

 // a map (int,int) -> int identifying the operator. use a flat_map for quicker access.
 std::map<std::pair<int, int>, int> int_pair_to_n;
 // boost::container::flat_map<std::pair<int, int>, int> int_pair_to_n;

};
}
