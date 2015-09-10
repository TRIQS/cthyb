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

namespace cthyb {

using namespace triqs::arrays;
using namespace triqs::gfs;
using namespace triqs::hilbert_space;
namespace h5 = triqs::h5;

using block_matrix_t = std::vector<matrix<double>>;
using full_hilbert_space_state_t = triqs::arrays::vector<double>;
using indices_t = fundamental_operator_set::indices_t;
using many_body_op_t = triqs::utility::many_body_operator<double>;
 
// Division of Hilbert Space into sub hilbert spaces, using the quantum numbers.
class sorted_spaces {

 using indices_t = fundamental_operator_set::indices_t;
 using many_body_op_t = triqs::operators::many_body_operator;
 using quantum_number_t = double;

 public:
 /// Not clean : why a vector, a std::vector ...
 struct eigensystem_t {
  vector<double> eigenvalues;    // in ascending order, the GS energy is set to 0 at initialisation
  std::vector<double> quantum_numbers;
  matrix<double> unitary_matrix; // H = U * \Lambda * U^+
 };

 sorted_spaces() = default;
 sorted_spaces(many_body_op_t const& h_, std::vector<many_body_op_t> const& qn_vector, fundamental_operator_set const& fops);
 sorted_spaces(many_body_op_t const& h_, fundamental_operator_set const& fops);

 sorted_spaces(many_body_op_t const& h_, std::vector<many_body_op_t> const& qn_vector)
    : sorted_spaces(h_, qn_vector, h_.make_fundamental_operator_set()) {}

 sorted_spaces(many_body_op_t const& h_) : sorted_spaces(h_, h_.make_fundamental_operator_set()) {}

 // useless ... fops is necessarly reordered !!!
 sorted_spaces(many_body_op_t const& h_, std::vector<many_body_op_t> const& qn_vector, std::vector<indices_t> const& indices_list)
    : sorted_spaces(h_, qn_vector, fops_from_indices(indices_list)) {}

 sorted_spaces(many_body_op_t const& h_, std::vector<indices_t> const& indices_list)
    : sorted_spaces(h_, fops_from_indices(indices_list)) {}

 /// Dimension of the full Hilbert space
 int dim() const { return _total_dim;}

 /// Number of Blocks
 int n_blocks() const { return eigensystems.size(); }

 /// The dimension of block b
 int get_block_dim(int b) const { return eigensystems[b].eigenvalues.size(); }

 /// Number of c operators
 int n_c_operators() const { return fops.size(); }

 /// FIX ME
 /// Connections for fundamental operators
 long fundamental_operator_connect(bool dagger, int index, int n) const {
  return (dagger ? creation_connection(index,n) : annihilation_connection(index,n));
 }

 /// FIX ME
 /// Connections for fundamental operators
 matrix<double> const& fundamental_operator_matrix(bool dagger, int index, int block_index) const {
  return (dagger ? cdag_matrices[index][block_index] : c_matrices[index][block_index]);
 }

 /// Ground state energy (i.e. min of all subspaces).
 double get_gs_energy() const { return gs_energy; }

 /// Get the i-th eigenvalue of block bl
 double get_eigenvalue(int bl, int i) const { return eigensystems[bl].eigenvalues[i]; }

 /// bl : block index, i : index within the block. Returns the index in the full hilbert space 
 int flatten_block_index(int bl, int i) const { return first_eigstate_of_block[bl] + i; }

 /// Vacuum is necessarly a block of size 1. Returns the block index.
 int get_vacuum_block_index() const { return _vacuum_index;}

 /// Vacuum is necessarly a block of size 1. Returns the block index.
 full_hilbert_space_state_t get_vacuum() const {
  full_hilbert_space_state_t st;
  st[flatten_block_index(_vacuum_index,0)] = 1;
  return st;
 }

 /// Trace (op * density_matrix)
 double average(block_matrix_t const& density_matrix, many_body_op_t const& op);

 /// Average projector
 double average_on_projector(block_matrix_t const& density_matrix, full_hilbert_space_state_t const& psi);

 /// Act with operator op on state st
 full_hilbert_space_state_t act(many_body_op_t const& op, full_hilbert_space_state_t const & st);

 /// Average the density matrix on the state psi

 /** 
  * For a symmetry implemented as a permutation of the C operators, return the matrix by block
  * Blocks must be invariant.
  * The permutation is a vector i-> j, where the int are the indices of the fundamental_operator_set of the object.
  */
 block_matrix_t matrix_element(std::vector<int> P);

 ///
 block_matrix_t matrix_element(std::vector<std::pair<indices_t, indices_t>> const &P);

 /// A human readable representation for eigenstate
 std::string eigenstate_repr(int bl, int k) const;

 ///
 std::vector<std::tuple<vector<double>,std::vector<double>>> get_energy_and_quantum_numbers() const;

 /// The partition function
 double partition_function(double beta) const;

 /// The atomic green function
 block_gf<imtime> atomic_gf(double beta, std::map<std::string, indices_t> const& indices_list, int n_tau,
                            std::vector<std::pair<int, int>> const& excluded_states = {}) const;

 /// The atomic observables
 //std::map<std::string,std::vector<double>> atomic_observables(std::map<std::string,many_body_op_t> const& obs_map) const;

 private:
 /// ------------------  DATA  -----------------

 // For each operator (index), the table of operator/connection to other sub_hilbert_spaces
 matrix<long> creation_connection, annihilation_connection;

 // For each operator (index), the matrix c_matrices[i][B] is the matrix from block B -> B'
 std::vector<std::vector<matrix<double>>> c_matrices, cdag_matrices;

 std::vector<sub_hilbert_space> sub_hilbert_spaces; // The subspace for each block, i.e. the list of fock states
 std::vector<eigensystem_t> eigensystems;           // Eigensystem in each subspace
 double gs_energy;                                  // Energy of the ground state
 fundamental_operator_set fops;                     // keep it to compute the local gf
 int _vacuum_index;                                 // Index of the bare vacuum

 // do not serialize. rebuild by complete_init2
 std::vector<int> first_eigstate_of_block; // Index of the first eigenstate of each block
 int _total_dim;
 std::vector<std::vector<quantum_number_t>> _quantum_numbers; // temporary variable used in construction phase. Do not USE.


 void autopartition(fundamental_operator_set const& fops, many_body_op_t const& h);
 void slice_hilbert_space_with_qn(many_body_op_t const& h_, std::vector<many_body_op_t> const& qn_vector,
                                  fundamental_operator_set const& fops);
 void complete_init1(many_body_op_t const& h_); // reorder the blocks, compute the matrices, ....
 void complete_init2();                         // final part also reused by h5_read

 // USED ONLY ONCE
 // Create matrix of an operator acting from one subspace to another (returns matrix + number of its nonzero elements)
 std::pair<matrix<double>, int> make_op_matrix(imperative_operator<hilbert_space> const& op, int from_sp, int to_sp) const;

 std::pair<int, matrix<double>> matrix_element_of_monomial(many_body_op_t::monomial_t const& op_vec, int B);
 friend std::ostream& operator<<(std::ostream& os, sorted_spaces const& ss);
 // h5
 friend std::string get_triqs_hdf5_data_scheme(sorted_spaces const&) { return "SortedSpaces"; }
 friend void h5_write(h5::group gr, std::string const& name, sorted_spaces const&);
 friend void h5_read(h5::group gr, std::string const& name, sorted_spaces&);

 // duplicate from solver_core : MOVE THIS
 fundamental_operator_set fops_from_indices(std::vector<indices_t> const& indices_list) {
  fundamental_operator_set fops;
  for (auto const& i : indices_list) {
   fops.insert(i);
  }
  return fops;
 }
};
}
