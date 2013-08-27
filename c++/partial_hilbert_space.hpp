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

#ifndef TRIQS_CTQMC_KRYLOV_HILBERT_SPACE_QN
#define TRIQS_CTQMC_KRYLOV_HILBERT_SPACE_QN

#include "fock_state.hpp"

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

class partial_hilbert_space {

  public:

  // constructor
  partial_hilbert_space() {}

  int n_bits (uint64_t f) const{ return (f * 0x200040008001ULL & 0x111111111111111ULL) % 0xf; }

  // add a fock state to the hilbert space basis
  void add_basis_fock(fock_state const & f) {
    size_t ind = fock_states.size();
    fock_states.push_back(f);
    fock_to_index.insert(std::make_pair(f, ind));
    fock_to_index_v.resize(f+1,-1);
    fock_to_index_v[f] = ind;

    //if (dimension() && (n_bits(*fock_states.begin()) != n_bits(f)))  TRIQS_RUNTIME_ERROR << "oops "<< f ;
  }

  // return dimension of the sub hilbert space
  size_t dimension() const { return fock_states.size(); }

  // find the index of a given state
  size_t get_state_index(fock_state f) const { 
   
   //std::cout << fock_to_index_v[f] << " "<< fock_to_index.find(f)->second<<std::endl;
   return fock_to_index_v[f];
   //return fock_to_index.find(f)->second; 
  }

  // return the state for a given index
  fock_state get_fock_state(size_t i) const { return fock_states[i]; }

  private:
      
  // the list of all fock states
  std::vector<fock_state> fock_states;

  // a map to quickly find the index of a state
  std::unordered_map<fock_state, size_t> fock_to_index;
  std::vector<long> fock_to_index_v;

};

}}}}
#endif
