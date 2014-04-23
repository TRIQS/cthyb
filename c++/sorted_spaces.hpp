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
#include "triqs/draft/hilbert_space_tools/state.hpp"

namespace cthyb {

using namespace triqs::arrays;
using namespace triqs::gfs;
using triqs::arrays::matrix;

// Division of Hilbert Space into sub hilbert spaces, using the quantum numbers.
class sorted_spaces {

 using indices_t = fundamental_operator_set::indices_t;
 using many_body_op_t = triqs::utility::many_body_operator<double>;

 public:
 struct eigensystem_t {
  vector<double> eigenvalues; // in ascending order, the GS energy is set to 0 at initialisation
  std::vector<state<sub_hilbert_space, double, false>> eigenstates;
  matrix<double> unitary_matrix; // H = U * \Lambda * U^+
 };

 // Constructor
 sorted_spaces() = default;
 sorted_spaces(many_body_op_t const& h_, std::vector<many_body_op_t> const& qn_vector, fundamental_operator_set const& fops);
 sorted_spaces(many_body_op_t const& h_, fundamental_operator_set const& fops);

 sorted_spaces(sorted_spaces const&) = delete;
 sorted_spaces(sorted_spaces&&) = default;
 sorted_spaces& operator=(sorted_spaces const&) = delete;
 sorted_spaces& operator=(sorted_spaces&&) = default;

 // Number of subspaces
 int n_subspaces() const { return sub_hilbert_spaces.size(); }

 // Number of c operators
 int n_c_operators() const { return fops.n_operators(); }

 // n-th subspace
 sub_hilbert_space const& subspace(int n) const { return sub_hilbert_spaces[n]; }

 // connections for fundamental operators
 long fundamental_operator_connect(bool dagger, int index, int n) const {
  return (dagger ? creation_connection[index][n] : annihilation_connection[index][n]);
 }

 // connections for fundamental operators
 matrix<double> const& fundamental_operator_matrix(bool dagger, int index, int block_index) const {
  return (dagger ? cdag_matrices[index][block_index] : c_matrices[index][block_index]);
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

 // For each operator (index), the table of operator/connection to other sub_hilbert_spaces
 std::vector<std::vector<long>> creation_connection, annihilation_connection;

 // For each operator (index), the matrice c_matrices[i][B] is the matrix from block B -> B'
 std::vector<std::vector<matrix<double>>> c_matrices, cdag_matrices;

 std::vector<eigensystem_t> eigensystems; // Eigensystem in each subspace
 double gs_energy;                        // Energy of the ground state
 fundamental_operator_set fops;           // keep it to compute the local gf

 void complete_init(many_body_op_t const& h_); // reorder the blocks, compute the matrices, ....
 void autopartition(fundamental_operator_set const& fops, many_body_op_t const& h);
 void slice_hilbert_space_with_qn(many_body_op_t const& h_, std::vector<many_body_op_t> const& qn_vector,
                                  fundamental_operator_set const& fops);
 friend std::ostream& operator<<(std::ostream& os, sorted_spaces const& ss);
};
}
