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
#pragma once
#include <map>
#include "state.hpp"
#include "sorted_spaces.hpp"
#include "hilbert_space.hpp"
#include <triqs/utility/time_pt.hpp>

namespace cthyb_krylov {

using triqs::utility::time_pt;

// The configuration of the Monte Carlo
struct configuration {

 struct op_desc {  // The description of the C operator
  int block_index; // the block index of the operator
  int inner_index; // the inner index inside the block
  bool dagger;     // is the operator a dagger
 };

 // a map associating an operator to an imaginary time
 using oplist_t=std::map<time_pt, op_desc, std::greater<time_pt>> ;
 // using oplist_t=boost::container::flat_map<time_pt, op_desc, std::greater<time_pt>> ;
 oplist_t oplist;

 // The boundary states, (subspace,state) pairs
 std::vector<std::pair<int, int>> boundary_block_states_ids;

 // construction and the basics stuff. value semantics, except = ?
 configuration(double beta_, sorted_spaces const& sosp, bool use_cutoff, double cutoff);

 double beta() const { return beta_; }

 friend std::ostream& operator<<(std::ostream& out, configuration const& c);

 template <class Archive> void serialize(Archive& ar, const unsigned int version) {
  ar& boost::serialization::make_nvp("oplist", oplist) & boost::serialization::make_nvp("beta", beta_) &
      boost::serialization::make_nvp("boundary_block_states_ids", boundary_block_states_ids);
 }

 private:
 double beta_;
};
}

