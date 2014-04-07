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

#include "./imperative_operator.hpp"
#include "./state.hpp"

using namespace triqs::arrays;

namespace cthyb_matrix {

 using namespace triqs::gfs;

// Division of Hilbert Space into sub hilbert spaces, using the quantum numbers.
class sorted_spaces {

 using indices_t = fundamental_operator_set::indices_t;

 public:
 struct eigensystem_t {
  vector<double> eigenvalues; // in ascending order, the GS energy is set to 0 at initialisation
  std::vector<state<sub_hilbert_space, double, false>> eigenstates;
  matrix<double> unitary_matrix; // H = U * \Lambda * U^+
 };

 // Constructor
 sorted_spaces(triqs::utility::many_body_operator<double> const& h_,
               std::vector<triqs::utility::many_body_operator<double>> const& qn_vector, fundamental_operator_set const& fops);

 sorted_spaces(sorted_spaces const&) = delete;
 sorted_spaces(sorted_spaces&&) = default;

 // Hamiltonian
 imperative_operator<sub_hilbert_space, false> const& get_hamiltonian() const { return hamiltonian; }

 // Number of subspaces
 int n_subspaces() const { return sub_hilbert_spaces.size(); }

 // Number of c operators
 int n_c_operators() const { return fops.n_operators(); }

 // n-th subspace
 sub_hilbert_space const& subspace(int n) const { return sub_hilbert_spaces[n]; }

 // an 0 state in n-th subpace
 state<sub_hilbert_space, double, false> substate(int n) const {
  return {subspace(n)};
 }

 // get fundamental operators
 //imperative_operator<sub_hilbert_space, true> const& get_fundamental_operator_from_linear_index(bool dagger,
 //                                                                                               int linear_index) const {
 // return (dagger ? creation_operators[linear_index] : destruction_operators[linear_index]);
// }

 // connections for fundamental operators
 long fundamental_operator_connect_from_linear_index(bool dagger, int linear_index, int n) const {
  return (dagger ? creation_connection[linear_index][n] : destruction_connection[linear_index][n]);
 }

 // connections for fundamental operators
 triqs::arrays::matrix<double> const& fundamental_operator_matrix_from_linear_index(bool dagger, int linear_index,
                                                                                    int block_index) const {
  return (dagger ? cdag_matrices[linear_index][block_index] : c_matrices[linear_index][block_index]);
 }

 // eigensystems for all blocks
 std::vector<eigensystem_t> const& get_eigensystems() const { return eigensystems; }

 // (global) gs energy (i.e. min of all subspaces).
 double get_gs_energy() const { return gs_energy; }

 /// The partition function
 double partition_function(double beta) const;

 /// The atomic green function
 block_gf<imtime> atomic_gf(double beta) const;


 private:
 /// ------------------  DATAS  -----------------

 std::vector<sub_hilbert_space> sub_hilbert_spaces; // all blocks

 imperative_operator<sub_hilbert_space, false> hamiltonian;

 // Given the linear index of the operator, return the table of operator/connection to other sub_hilbert_spaces
 std::vector<std::vector<long>> creation_connection, destruction_connection;
 //std::vector<imperative_operator<sub_hilbert_space, true>> creation_operators, destruction_operators;

 // Given the linear index of the operator i, the matrice c_matrices[i][B] is the matrix from block B -> B'
 std::vector<std::vector<triqs::arrays::matrix<double>>> c_matrices, cdag_matrices;

 // Eigensystem in each subspace
 std::vector<eigensystem_t> eigensystems;

 // Energy of the ground state
 double gs_energy;

 // keep it to compute the local gf
 fundamental_operator_set fops;

 void complete_init(triqs::utility::many_body_operator<double> const& h_); // reorder the blocks, compute the matrices, ....

 void slice_hilbert_space_with_qn(triqs::utility::many_body_operator<double> const& h_,
                                 std::vector<triqs::utility::many_body_operator<double>> const& qn_vector,
                                 fundamental_operator_set const& fops);

 friend std::ostream& operator<<(std::ostream& os, sorted_spaces const& ss);
};
}
