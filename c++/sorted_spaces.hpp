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
#include <triqs/arrays.hpp>
#include <string>
#include <vector>
#include <map>
//#include <triqs/utility/tuple_tools.hpp>
#include <triqs/arrays/linalg/eigenelements.hpp>

#include "fundamental_operator_set.hpp"
#include "hilbert_space.hpp"
#include "operator.hpp"
#include "imperative_operator.hpp"
#include "state.hpp"

using namespace triqs::arrays;

namespace cthyb_krylov {

/// ???
struct block_desc_t {
 std::string name;
 std::vector<fundamental_operator_set::indices_t> indices;
 // for python interface only
 void indices_push_back(std::string a, std::string b) {
  indices.push_back({a, b});
 }
};

// Division of Hilbert Space into sub hilbert spaces, using the quantum numbers.
class sorted_spaces {

 using indices_t = fundamental_operator_set::indices_t;
 using quantum_number_t = double;
 
 public:
 struct eigensystem_t {
  vector<double> eigenvalues; // in ascending order, the GS energy is set to 0 at initialisation
  std::vector<state<sub_hilbert_space, double, false>> eigenstates;
  matrix<double> unitary_matrix; // H = U * \Lambda * U^+
 };

 // Constructor
 sorted_spaces(triqs::utility::many_body_operator<double> const& h_,
               std::vector<triqs::utility::many_body_operator<double>> const& qn_vector, fundamental_operator_set const& fops,
               std::vector<block_desc_t> const& block_structure);

 // Hamiltonian
 imperative_operator<sub_hilbert_space, false> const& get_hamiltonian() const { return hamilt; }

 // Number of subspaces
 int n_subspaces() const { return n_blocks; }

 // n-th subspace
 sub_hilbert_space const& subspace(int n) const { return *sub_hilbert_spaces[n]; }

 // an 0 state in n-th subpace
 state<sub_hilbert_space, double, false> substate(int n) const {
  return {subspace(n)};
 }

 // get fundamental operators
 imperative_operator<sub_hilbert_space, true> const& get_fundamental_operator(bool dagger, int block_index,
                                                                              int inner_index) const {
  auto p = int_pair_to_n.at({block_index, inner_index});
  return (dagger ? creation_operators[p] : destruction_operators[p]);
 }

 // connections for fundamental operators
 long fundamental_operator_connect(bool dagger, int block_index, int inner_index, int n) const {
  auto p = int_pair_to_n.at({block_index, inner_index});
  return (dagger ? creation_connection[p][n] : destruction_connection[p][n]);
 }

 // eigensystems for all blocks
 std::vector<eigensystem_t> const& get_eigensystems() const { return eigensystems; }

 // (global) gs energy (i.e. min of all subspaces).
 double get_gs_energy() const { return gs_energy; }

 private:

 /// ------------------  DATAS  -----------------
 
 int n_blocks; // number of sub_hilbert_spaces

 // a map (int,int) -> int identifying the operator
 std::map<std::pair<int, int>, int> int_pair_to_n;

 // the hamiltonian,
 imperative_operator<sub_hilbert_space, false> hamilt;

// have to use shared_ptr, because sub_hilbert_space is non-copyable
 std::vector<std::shared_ptr<sub_hilbert_space>> sub_hilbert_spaces;
 //std::vector<sub_hilbert_space> sub_hilbert_spaces;

 // Given the index of a ???? EXPLAiN
 std::vector<std::vector<long>> creation_connection, destruction_connection;
 std::vector<imperative_operator<sub_hilbert_space, true>> creation_operators, destruction_operators;

 // Eigensystem in each subspace
 std::vector<eigensystem_t> eigensystems;

 // Energy of the ground state
 double gs_energy;

 // Keep the QN, only for later printing ? OR MAKE THE STRING ...
 std::vector<std::vector<quantum_number_t>> quantum_numbers;
 
 /// ------------------  Functions  -----------------
 
 friend std::ostream& operator<<(std::ostream& os, sorted_spaces const& ss);

 void compute_eigensystems(); // auxiliary function to compute the eigensystem
};
}
