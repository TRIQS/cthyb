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

#include "fock_state.hpp"
#include <boost/container/flat_map.hpp>

namespace cthyb_krylov {

// a subhilbert space, as a set of basis Fock states.
// contains 2 functions to switch from the Fock state to its number in this set
class partial_hilbert_space {

 public:
 partial_hilbert_space(int index) : index(index) {}

 // add a fock state to the hilbert space basis
 void add_fock_state(fock_state f) {
  int ind = fock_states.size();
  fock_states.push_back(f);
  fock_to_index.insert(std::make_pair(f, ind));
 }

 // dimension of the sub hilbert space
 int dimension() const { return fock_states.size(); }

 // find the index of a given state
 int get_state_index(fock_state f) const { return fock_to_index.find(f)->second; }

 // the state for a given index
 fock_state get_fock_state(int i) const { return fock_states[i]; }

 int get_index() const {
  return index;
 };

 private:
 // index of a partial space
 int index;

 // the list of all fock states
 std::vector<fock_state> fock_states;

 // reverse a map to quickly find the index of a state
 // the boost flat_map is implemented as an ordered vector,
 // hence is it slow to insert (we don't care) but fast to look up (we do it a lot)
 std::map<fock_state, int> fock_to_index;
 //boost::container::flat_map<fock_state, int> fock_to_index;
};
}
