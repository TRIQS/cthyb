#pragma once
#include <triqs/arrays.hpp>
#include <triqs/gfs.hpp>
#include <string>
#include <vector>
#include <map>
#include <triqs/arrays/linalg/eigenelements.hpp>

#include "./imperative_operator.hpp"
#include "./state.hpp"

using namespace triqs::arrays;

namespace cthyb_krylov {

 using namespace triqs::gfs;

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
 long get_fundamental_operator_linear_index(int block_index, int inner_index) const {
  return int_pair_to_n.at({block_index, inner_index});
 }

 // get fundamental operators
 imperative_operator<sub_hilbert_space, true> const& get_fundamental_operator(bool dagger, int block_index,
                                                                              int inner_index) const {
  auto p = int_pair_to_n.at({block_index, inner_index});
  return (dagger ? creation_operators[p] : destruction_operators[p]);
 }

 // get fundamental operators
 imperative_operator<sub_hilbert_space, true> const& get_fundamental_operator_from_linear_index(bool dagger,
                                                                                                int linear_index) const {
  return (dagger ? creation_operators[linear_index] : destruction_operators[linear_index]);
 }

 // connections for fundamental operators
 long fundamental_operator_connect(bool dagger, int block_index, int inner_index, int n) const {
  auto p = int_pair_to_n.at({block_index, inner_index});
  return (dagger ? creation_connection[p][n] : destruction_connection[p][n]);
 }

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

 /// The atomic green function
 block_gf<imtime> atomic_gf(double beta) const;

 private:
 /// ------------------  DATAS  -----------------

 std::vector<sub_hilbert_space> sub_hilbert_spaces; // all blocks

 // a map (int,int) -> int identifying the operator. use a flat_map for quicker access.
 std::map<std::pair<int, int>, int> int_pair_to_n;
 // boost::container::flat_map<std::pair<int, int>, int> int_pair_to_n;

 imperative_operator<sub_hilbert_space, false> hamiltonian;

 // Given the linear index of the operator, return the table of operator/connection to other sub_hilbert_spaces
 std::vector<std::vector<long>> creation_connection, destruction_connection;
 std::vector<imperative_operator<sub_hilbert_space, true>> creation_operators, destruction_operators;

 // Given the linear index of the operator i, the matrice c_matrices[i][B] is the matrix from block B -> B'
 std::vector<std::vector<triqs::arrays::matrix<double>>> c_matrices, cdag_matrices;

 // Eigensystem in each subspace
 std::vector<eigensystem_t> eigensystems;

 // Energy of the ground state
 double gs_energy;

 // Keep the QN, only for later printing ? OR MAKE THE STRING ...
 std::vector<std::vector<quantum_number_t>> quantum_numbers;

 // keep it to compute the local gf
 fundamental_operator_set fops;

 friend std::ostream& operator<<(std::ostream& os, sorted_spaces const& ss);
};
}
