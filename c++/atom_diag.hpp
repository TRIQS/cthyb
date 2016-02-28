/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/state.hpp>
#include <triqs/hilbert_space/imperative_operator.hpp>
#include "./array_suppl.hpp"

namespace cthyb {

using namespace triqs::arrays;
using namespace triqs::gfs;
using namespace triqs::hilbert_space;
namespace h5 = triqs::h5;
namespace operators = triqs::operators;

using h_scalar_t = double;                                               // type of scalar for H_loc: double or complex.
using block_matrix_t = std::vector<matrix<h_scalar_t>>;                  // block diagonal matrix
using full_hilbert_space_state_t = triqs::arrays::vector<h_scalar_t>;    // the big vector in the full Hilbert space
using indices_t = fundamental_operator_set::indices_t;                   //
using many_body_op_t = triqs::operators::many_body_operator_real;
using quantum_number_t = double;

// Division of Hilbert Space into sub hilbert spaces, using the quantum numbers.
class atom_diag {
 friend class atom_diag_worker;

 public:
 struct eigensystem_t {
  triqs::arrays::vector<double> eigenvalues; // in ascending order, the GS energy is set to 0 at initialisation
  matrix<h_scalar_t> unitary_matrix;         // H = U * \Lambda * U^+, from the Fock space basis to the block basis
 };

 TRIQS_CPP2PY_IGNORE atom_diag() = default;
 atom_diag(many_body_op_t const& h_, fundamental_operator_set const& fops);
 atom_diag(many_body_op_t const& h_, fundamental_operator_set const& fops, std::vector<many_body_op_t> const& qn_vector);

 /// The Hamiltonian
 many_body_op_t const& get_h_atomic() const { return h_atomic; }

 /// The fundamental operator set used at construction
 fundamental_operator_set const& get_fops() const { return fops; }

 /// Dimension of the full Hilbert space
 int get_full_hilbert_space_dim() const { return _total_dim; }

 /// Number of Blocks
 int n_blocks() const { return eigensystems.size(); }

 /// The dimension of block b
 int get_block_dim(int b) const { return eigensystems[b].eigenvalues.size(); }

 /// Returns the index in the full hilbert space for block_index and i, the index within the block.
 int flatten_block_index(int block_index, int i) const { return first_eigstate_of_block[block_index] + i; }

 /// Return the range of indices of block B
 TRIQS_CPP2PY_IGNORE range index_range_of_block(int bl) const {
  return range{first_eigstate_of_block[bl], first_eigstate_of_block[bl] + get_block_dim(bl) + 1};
 }

 /// Get the eigensystem
 TRIQS_CPP2PY_IGNORE std::vector<eigensystem_t> const& get_eigensystem() const { return eigensystems; }

 /// Get the i-th eigenvalue of block bl
 double get_eigenvalue(int block_index, int i) const { return eigensystems[block_index].eigenvalues[i]; }

 /// A vector of all the energies, by blocks. result[block_number][i] is the energy
 std::vector<std::vector<double>> get_energies() const;

 /// A vector of all the QNs, by blocks : result[block_number][qn_index] is the .....
 std::vector<std::vector<double>> const& get_quantum_numbers() const { return quantum_numbers; }

 /// Ground state energy (i.e. min of all subspaces).
 double get_gs_energy() const { return gs_energy; }

 /// Returns the block index of the vacuum state.
 int get_vacuum_block_index() const { return vacuum_block_index; }

 /// Returns the inner index of the vacuum state.
 int get_vacuum_inner_index() const { return vacuum_inner_index; }

 /// Returns the vacuum state as a long vector in the full Hilbert space.
 full_hilbert_space_state_t get_vacuum_state() const;

 /**
 * Connections for fundamental operators C
 *
 * op_linear_index : the linear index (i.e. number) of the c operator, as defined by the fundamental_operator_set fops
 * block_number : the number of the initial block
 * @return : the number of the final block
 */
 long c_connection(int op_linear_index, int block_index) const { return annihilation_connection(op_linear_index, block_index); }
 /**
   * Connections for fundamental operators C^\dagger
   *
   * op_linear_index : the linear index (i.e. number) of the c operator, as defined by the fundamental_operator_set fops
   * block_number : the number of the initial block
   * @return : the number of the final block
   */
 long cdag_connection(int op_linear_index, int block_index) const { return creation_connection(op_linear_index, block_index); }

 /**
  * Matrix for fundamental operators C
  *
  * op_linear_index : the linear index (i.e. number) of the c operator, as defined by the fundamental_operator_set fops
  * block_number : the number of the initial block
  * @return : the number of the final block
  */
 matrix<h_scalar_t> const& c_matrix(int op_linear_index, int block_index) const {
  return c_matrices[op_linear_index][block_index];
 }

 /**
  * Matrix for fundamental operators C^\dagger
  *
  * op_linear_index : the linear index (i.e. number) of the c operator, as defined by the fundamental_operator_set fops
  * block_number : the number of the initial block
  * @return : the number of the final block
  */
 matrix<h_scalar_t> const& cdag_matrix(int op_linear_index, int block_index) const {
  return cdag_matrices[op_linear_index][block_index];
 }

 //FIXME
 /**
   * For a symmetry S implemented as a unitary transformation of the C operators
   * returns the matrix representation of S, which is block diagonal.
   * Blocks must be invariant by the symmetry of the function throws.
   */
 // block_matrix_t matrix_of_symmetry(matrix<h_scalar_t> const& S);

 /**
  * Given a monomial (ccccc), and a block B, returns
  *  - the block connected by ccccc from B
  *  - the corresponding block matrix (not necessarly square)
  */
 TRIQS_CPP2PY_IGNORE std::pair<int, matrix<h_scalar_t>> matrix_element_of_monomial(operators::monomial_t const& op_vec, int B) const;

 private:
 /// ------------------  DATA  -----------------

 many_body_op_t h_atomic;                                    // Atomic hamiltonian
 fundamental_operator_set fops;                              // keep it to compute the local gf
 std::vector<sub_hilbert_space> sub_hilbert_spaces;          // The subspace for each block, i.e. the list of fock states
 std::vector<eigensystem_t> eigensystems;                    // Eigensystem in each subspace
 matrix<long> creation_connection, annihilation_connection;  // creation_connection[operator_linear_index][B] -> B', final block
 std::vector<std::vector<matrix<h_scalar_t>>> c_matrices;    // c_matrices[operator_linear_index][B] = matrix from block B -> B'
 std::vector<std::vector<matrix<h_scalar_t>>> cdag_matrices; // idem for c dagger operators
 double gs_energy;                                           // Energy of the ground state
 int vacuum_block_index;                                     // Block index of the bare vacuum
 int vacuum_inner_index;                                     // Inner index of the bare vacuum

 std::vector<std::vector<h_scalar_t>> quantum_numbers; // values of the quantum number for this blocks

 // do not serialize. rebuild by complete_init
 void complete_init();
 std::vector<int> first_eigstate_of_block; // Index of the first eigenstate of each block
 int _total_dim;                           // total_dimension of the Hilbert_space

 friend std::ostream& operator<<(std::ostream& os, atom_diag const& ss);
 friend std::string get_triqs_hdf5_data_scheme(atom_diag const&) { return "AtomicDiagonalization"; }
 friend void h5_write(h5::group gr, std::string const& name, atom_diag const&);
 friend void h5_read(h5::group gr, std::string const& name, atom_diag&);
};
}
