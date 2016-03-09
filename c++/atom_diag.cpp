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
#include "./atom_diag.hpp"
#include "./atom_diag_worker.hpp"
#include <sstream>
#include <bitset>
#include <algorithm>
using namespace triqs::arrays;
using namespace triqs::hilbert_space;
using std::string;

namespace cthyb {

atom_diag::atom_diag(many_body_op_t const& h_, fundamental_operator_set const& fops, std::vector<many_body_op_t> const& qn_vector)
   : h_atomic(h_), fops(fops) {
 atom_diag_worker{this}.partition_with_qn(qn_vector);
 complete_init();
}

//-----------------------------

atom_diag::atom_diag(many_body_op_t const& h_, fundamental_operator_set const& fops) : h_atomic(h_), fops(fops) {
 atom_diag_worker{this}.autopartition();
 complete_init();
}

// -----------------------------------------------------------------

void atom_diag::complete_init() {
 _total_dim = 0;
 for (auto const& es : eigensystems) _total_dim += es.eigenvalues.size();

 // Calculate the index of the first eigenstate of each block
 first_eigstate_of_block.resize(_total_dim, 0);
 for (int bl = 1; bl < n_blocks(); ++bl) first_eigstate_of_block[bl] = first_eigstate_of_block[bl - 1] + get_block_dim(bl - 1);
}

// -----------------------------------------------------------------
std::vector<std::vector<double>> atom_diag::get_energies() const {
 std::vector<std::vector<double>> R;
 for (auto const& es : eigensystems) { 
  std::vector<double> v(es.eigenvalues.size());
  for (int i = 0; i < es.eigenvalues.size(); ++i) v[i] = es.eigenvalues[i];
  R.push_back(v);
 }
 return R;
}

// -----------------------------------------------------------------

full_hilbert_space_state_t atom_diag::get_vacuum_state() const {
 full_hilbert_space_state_t st(_total_dim);
 st() = 0;
 st[flatten_block_index(vacuum_block_index, vacuum_inner_index)] = 1;
 return st;
}

// -----------------------------------------------------------------

std::ostream& operator<<(std::ostream& os, atom_diag const& ss) {

 os << "Dimension of full Hilbert space: " << ss.get_full_hilbert_space_dim() << std::endl;
 os << "Number of blocks: " << ss.n_blocks() << std::endl;
 for (int n = 0; n < ss.n_blocks(); ++n) {
  os << "Block " << n << ", ";
  os << "relative gs energy : " << ss.eigensystems[n].eigenvalues[0] << std::endl;
  os << "size = " << ss.eigensystems[n].eigenvalues.size() << std::endl;
  os << "-------------------------" << std::endl;
 }
 return os;
}

// -----------------------------------------------------------------

namespace {

 // We assume i<j. S = change of sign
 std::pair<bool, fock_state_t> act_with_transposition(fock_state_t f, int i, int j, bool S) {
  std::bitset<32> fb(f);
  bool i_set = fb[i], j_set = fb[j];
  if (!i_set && !j_set) return {S, f}; // no bit sets, no effect
  if (i_set && j_set) return {!S, f};  // both bit sets, need to exchange the c^dagger.
  for (int u = i + 1; u < j; ++u)
   if (fb[u]) S = !S;
  if (!i_set) std::swap(i, j);
  fb[i] = 0;
  fb[j] = 1;
  return {S, fb.to_ullong()};
 }
}

//---------------------
/*
FIXME
block_matrix_t atom_diag::matrix_element(std::vector<int> P) {

 int L = P.size();
 TRIQS_ASSERT(L == fops.size()) << "The size of the permutation is incorrect " << L << "  " << fops.size();

 // Decompose the permutation into transposition
 std::vector<std::pair<int, int>> transpositions;
 for (int i = 0; i < L; ++i) {
  int k = P[i];
  if (k == i) continue;
  // we use the transposition i <-> k
  transpositions.push_back({i, k});
  // product T_ik * P : i -> k,  k -> i
  for (int j = i + 1; j < L; ++j)
   if (P[j] == i) {
    P[j] = k;
    break;
   }
 }

 // For each block, we build the matrix of the permutation
 block_matrix_t res;

 for (int bl = 0; bl < n_blocks(); ++bl) {
  matrix<h_scalar_t> R(get_block_dim(bl), get_block_dim(bl));
  R() = 0;
  auto const& SP = sub_hilbert_spaces[bl];
  for (auto fock : SP.get_all_fock_states()) {
   bool S = false;
   fock_state_t f = fock;
   for (auto const& T : transpositions) {
    std::tie(S, f) = act_with_transposition(f, T.first, T.second, S);
   }
   TRIQS_ASSERT(SP.has_state(f)) << " State f " << std::bitset<32>(f) << " is not in the sub_hilbert_space";
   R(SP.get_state_index(f), SP.get_state_index(fock)) = (S ? -1 : 1);
  }

  // Now we go to the eigenstate basis
  auto const& U = eigensystems[bl].unitary_matrix;
  R = U.transpose() * R * U; // COMPLEX CASE : change to dagger !
  for (auto& x : R)          // clean for clear printing
   if (std::abs(x) < 1.e-10) x = 0;
  res.push_back(R);
 }
 return res;
}

//---------------------

block_matrix_t atom_diag::matrix_element(std::vector<std::pair<indices_t, indices_t>> const& P) {
 int L = fops.size();
 std::vector<int> Pi(L);
 for (int i = 0; i < L; ++i) Pi[i] = i;
 for (auto const& ii : P) Pi[fops[ii.first]] = fops[ii.second];
 return matrix_element(Pi);
}
*/
// -----------------------------------------------------------------
// FIXME MOVE OUT OF CLASS
std::pair<int, matrix<h_scalar_t>> atom_diag::matrix_element_of_monomial(operators::monomial_t const& op_vec, int B) const {

 matrix<h_scalar_t> m = triqs::arrays::make_unit_matrix<h_scalar_t>(get_block_dim(B));
 for (int i = op_vec.size() - 1; i >= 0; --i) {
  int ind = fops[op_vec[i].indices];
  int Bp = (op_vec[i].dagger ? creation_connection(ind, B) : annihilation_connection(ind, B));
  if (Bp == -1) return {-1, std::move(m)};
  m = (op_vec[i].dagger ? cdag_matrices[ind][B] : c_matrices[ind][B]) * m;
  B = Bp;
 }
 return {B, std::move(m)};
}

}

// FIXME move into the library
namespace triqs {
namespace hilbert_space {
 // -----------------------------------------------------------------
 std::string get_triqs_hdf5_data_scheme(sub_hilbert_space const&) { return "sub_hilbert_space"; }

 // -----------------------------------------------------------------

 void h5_write(h5::group fg, std::string const& name, sub_hilbert_space const& x) {
  auto gr = fg.create_group(name);
  h5_write(gr, "fock_states", x.get_all_fock_states());
  h5_write(gr, "index", x.get_index());
 }

 // -----------------------------------------------------------------
 void h5_read(h5::group fg, std::string const& name, sub_hilbert_space& x) {
  using h5::h5_read;
  auto gr = fg.open_group(name);
  auto fs = h5_read<std::vector<fock_state_t>>(gr, "fock_states");
  auto index = h5_read<int>(gr, "index");
  x = sub_hilbert_space{index};
  for (auto const& s : fs) x.add_fock_state(s);
 }
}
}

namespace cthyb {
// -----------------------------------------------------------------
std::string get_triqs_hdf5_data_scheme(atom_diag::eigensystem_t const&) { return "atom_diag::eigensystem_t"; }

// -----------------------------------------------------------------

void h5_write(h5::group fg, std::string const& name, atom_diag::eigensystem_t const& x) {
 auto gr = fg.create_group(name);
 h5_write(gr, "eigenvalues", x.eigenvalues);
 h5_write(gr, "unitary_matrix", x.unitary_matrix);
}

// -----------------------------------------------------------------
void h5_read(h5::group fg, std::string const& name, atom_diag::eigensystem_t& x) {
 auto gr = fg.open_group(name);
 h5_read(gr, "eigenvalues", x.eigenvalues);
 h5_read(gr, "unitary_matrix", x.unitary_matrix);
}

// -----------------------------------------------------------------
void h5_write(h5::group fg, std::string const& name, atom_diag const& x) {
 auto gr = fg.create_group(name);
 gr.write_triqs_hdf5_data_scheme(x);
 h5_write(gr, "creation_connection", x.creation_connection);
 h5_write(gr, "annihilation_connection", x.annihilation_connection);

 auto write_sparse = [&](std::string na, std::vector<std::vector<matrix<h_scalar_t>>> const& Mvv) {
  auto gr2 = gr.create_group(na);
  for (int i = 0; i < Mvv.size(); ++i)
   for (int j = 0; j < Mvv[i].size(); ++j)
    if (!Mvv[i][j].is_empty()) h5_write(gr2, std::to_string(i) + ' ' + std::to_string(j), Mvv[i][j]);
 };

 write_sparse("c_matrices", x.c_matrices);
 write_sparse("cdag_matrices", x.cdag_matrices);

 h5_write(gr, "sub_hilbert_spaces", x.sub_hilbert_spaces);
 h5_write(gr, "eigensystems", x.eigensystems);
 h5_write(gr, "gs_energy", x.gs_energy);
 h5_write_attribute(gr, "fops", x.fops);
 h5_write(gr, "vacuum_block_index", x.vacuum_block_index);
 h5_write(gr, "vacuum_inner_index", x.vacuum_inner_index);
}

// -----------------------------------------------------------------
void h5_read(h5::group fg, std::string const& name, atom_diag& x) {
 auto gr = fg.open_group(name);
 h5_read(gr, "creation_connection", x.creation_connection);
 h5_read(gr, "annihilation_connection", x.annihilation_connection);

 auto read_sparse = [&](std::string na, std::vector<std::vector<matrix<h_scalar_t>>>& Mvv) {
  Mvv.resize(first_dim(x.creation_connection), std::vector<matrix<h_scalar_t>>(second_dim(x.creation_connection)));
  auto gr2 = gr.open_group(na);
  for (auto s : gr2.get_all_dataset_names()) {
   std::stringstream ss(s);
   std::string item1, item2;
   std::getline(ss, item1, ' ');
   std::getline(ss, item2, ' ');
   int i = std::stoi(item1), j = std::stoi(item2);
   h5_read(gr2, s, Mvv[i][j]);
  }
 };

 read_sparse("c_matrices", x.c_matrices);
 read_sparse("cdag_matrices", x.cdag_matrices);

 h5_read(gr, "sub_hilbert_spaces", x.sub_hilbert_spaces);
 h5_read(gr, "eigensystems", x.eigensystems);
 h5_read(gr, "gs_energy", x.gs_energy);
 h5_read_attribute(gr, "fops", x.fops);
 h5_read(gr, "vacuum_block_index", x.vacuum_block_index);
 h5_read(gr, "vacuum_inner_index", x.vacuum_inner_index);
 x.complete_init();
}
}
