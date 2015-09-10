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
#include "./sorted_spaces.hpp"
#include <triqs/arrays/linalg/eigenelements.hpp>
#include <triqs/hilbert_space/space_partition.hpp>
#include <sstream>
#include <bitset>
#include <algorithm>
using namespace triqs::arrays;
using namespace triqs::hilbert_space;
using std::string;

namespace cthyb {

/// CHANGE FOR ONE CONSTRUCTOR ONLY + string for option
sorted_spaces::sorted_spaces(many_body_op_t const& h_, std::vector<many_body_op_t> const& qn_vector,
                             fundamental_operator_set const& fops)
   : fops(fops) {
 slice_hilbert_space_with_qn(h_, qn_vector, fops);
 complete_init1(h_);
}

//-----------------------------

sorted_spaces::sorted_spaces(many_body_op_t const& h_, fundamental_operator_set const& fops) : fops(fops) {
 autopartition(fops, h_);
 complete_init1(h_);
}

//---------------------------------------------------------------------------------

void sorted_spaces::autopartition(fundamental_operator_set const& fops, many_body_op_t const& h) {

 hilbert_space full_hs(fops);

 imperative_operator<hilbert_space, double, false> hamiltonian(h, fops);
 state<hilbert_space, double, true> st(full_hs);

 using space_partition_t = space_partition<state<hilbert_space, double, true>, imperative_operator<hilbert_space, double, false>>;
 // Split the Hilbert space
 space_partition_t SP(st, hamiltonian, false);

 std::vector<space_partition_t::matrix_element_map_t> creation_melem(fops.size());
 std::vector<space_partition_t::matrix_element_map_t> annihilation_melem(fops.size());
 // Merge subspaces
 for (auto const& o : fops) {
  auto create = many_body_op_t::make_canonical(true, o.index);
  auto destroy = many_body_op_t::make_canonical(false, o.index);

  imperative_operator<hilbert_space, double> op_c_dag(create, fops), op_c(destroy, fops);

  int n = o.linear_index;
  std::tie(creation_melem[n], annihilation_melem[n]) = SP.merge_subspaces(op_c_dag, op_c, true);
 }

 // Fill subspaces
 sub_hilbert_spaces.reserve(SP.n_subspaces());
 for (int n = 0; n < SP.n_subspaces(); ++n) sub_hilbert_spaces.emplace_back(n);

 foreach(SP, [&](fock_state_t s, int spn) { sub_hilbert_spaces[spn].add_fock_state(s); });

 // Fill connections
 creation_connection.resize(fops.size(), SP.n_subspaces());
 annihilation_connection.resize(fops.size(), SP.n_subspaces());
 creation_connection.as_array_view() = -1;
 annihilation_connection.as_array_view() = -1;

 for (auto const& o : fops) {
  int n = o.linear_index;
  for (auto const& e : creation_melem[n]) {
   fock_state_t i, f;
   std::tie(i, f) = e.first;
   creation_connection(n,SP.lookup_basis_state(i)) = SP.lookup_basis_state(f);
  }
  for (auto const& e : annihilation_melem[n]) {
   fock_state_t i, f;
   std::tie(i, f) = e.first;
   annihilation_connection(n,SP.lookup_basis_state(i)) = SP.lookup_basis_state(f);
  }
 }
}

//-----------------------------

// define a more tolerant comparison between vectors for the quantum numbers
struct lt_dbl {
 bool operator()(std::vector<double> const& v1, std::vector<double> const& v2) const {
  for (int i = 0; i < v1.size(); ++i) {
   if (v1[i] < (v2[i] - 1e-8))
    return true;
   else if (v2[i] < (v1[i] - 1e-8))
    return false;
  }
  return false;
 }
};

void sorted_spaces::slice_hilbert_space_with_qn(many_body_op_t const& h_, std::vector<many_body_op_t> const& qn_vector,
                                                fundamental_operator_set const& fops) {

 hilbert_space full_hs(fops);

 // hilbert spaces and quantum numbers
 std::map<std::vector<double>, int, lt_dbl> map_qn_n;

 // The QN as operators : a vector of imperative operators for the quantum numbers
 std::vector<imperative_operator<hilbert_space, double>> qsize;
 for (auto& qn : qn_vector) qsize.emplace_back(qn, fops);

 // Helper function to get quantum numbers
 auto get_quantum_numbers = [&qsize](state<hilbert_space, double, false> const& s) {
  std::vector<quantum_number_t> qn;
  for (auto const& op : qsize) qn.push_back(dot_product(s, op(s)));
  return qn;
 };

 /*
   The first part consists in dividing the full Hilbert space
   into smaller subspaces using the quantum numbers
 */
 for (int r = 0; r < full_hs.size(); ++r) {

  // fock_state corresponding to r
  fock_state_t fs = full_hs.get_fock_state(r);

  // the state we'll act on
  state<hilbert_space, double, false> s(full_hs);
  s(r) = 1.0;

  // create the vector with the quantum numbers
  std::vector<quantum_number_t> qn = get_quantum_numbers(s);

  // if first time we meet these quantum numbers create partial Hilbert space
  if (map_qn_n.count(qn) == 0) {
   auto n_blocks = sub_hilbert_spaces.size();
   sub_hilbert_spaces.emplace_back(n_blocks); // a new sub_hilbert_space
   _quantum_numbers.push_back(qn);
   map_qn_n[qn] = n_blocks;
  }

  // add fock state to partial Hilbert space
  sub_hilbert_spaces[map_qn_n[qn]].add_fock_state(fs);
 }

 // ----  now make the creation map -----

 auto creation_map = std::vector<std::vector<int>>(fops.size(), std::vector<int>(sub_hilbert_spaces.size(), -1));
 auto annihilation_map = creation_map;

 // init the mapping tables
 creation_connection.resize(fops.size(), sub_hilbert_spaces.size());
 annihilation_connection.resize(fops.size(), sub_hilbert_spaces.size());
 creation_connection.as_array_view() = -1;
 annihilation_connection.as_array_view() = -1;

 for (auto const& x : fops) {

  // get the operators and their index
  int n = x.linear_index;
  auto create = many_body_op_t::make_canonical(true, x.index);
  auto destroy = many_body_op_t::make_canonical(false, x.index);

  // construct their imperative counterpart
  imperative_operator<hilbert_space, double> op_c_dag(create, fops), op_c(destroy, fops);

  // to avoid declaring every time in the loop below
  std::vector<quantum_number_t> qn_before, qn_after;

  // now act on the state with the c, c_dag to see how quantum numbers change
  for (int r = 0; r < full_hs.size(); ++r) {

   // the state we'll act on and its quantum numbers
   state<hilbert_space, double, false> s(full_hs);
   s(r) = 1.0;
   qn_before = get_quantum_numbers(s);

   // apply creation on state to figure quantum numbers
   qn_after = get_quantum_numbers(op_c_dag(s));

   // insert in creation map checking whether it was already there
   if (!triqs::utility::is_zero(dot_product(op_c_dag(s), op_c_dag(s)))) {
    auto origin = sub_hilbert_spaces[map_qn_n[qn_before]].get_index();
    auto target = sub_hilbert_spaces[map_qn_n[qn_after]].get_index();
    if (creation_map[n][origin] == -1)
     creation_map[n][origin] = target;
    else if (creation_map[n][origin] != target)
     TRIQS_RUNTIME_ERROR << "Internal Error, Sorted Space, creation";
    creation_connection(n,map_qn_n[qn_before]) = map_qn_n[qn_after];
   }

   // apply annihilation on state to figure quantum numbers
   qn_after = get_quantum_numbers(op_c(s));

   // insert in annihilation map checking whether it was already there
   if (!triqs::utility::is_zero(dot_product(op_c(s), op_c(s)))) {
    auto origin = sub_hilbert_spaces[map_qn_n[qn_before]].get_index();
    auto target = sub_hilbert_spaces[map_qn_n[qn_after]].get_index();
    if (annihilation_map[n][origin] == -1)
     annihilation_map[n][origin] = target;
    else if (annihilation_map[n][origin] != target)
     TRIQS_RUNTIME_ERROR << "Internal Error, Sorted Space, annihilation";
    annihilation_connection(n,map_qn_n[qn_before]) = map_qn_n[qn_after];
   }
  }
 }
}

// -------------------------------------------------------------------------------------------------

//// WHY DO WE RETURN NON_ZEROS : useless

std::pair<matrix<double>, int> sorted_spaces::make_op_matrix(imperative_operator<hilbert_space> const& op, int from_spn,
                                                             int to_spn) const {

 hilbert_space full_hs(fops);
 int nonzeros = 0;
 auto const& from_sp = sub_hilbert_spaces[from_spn];
 auto const& to_sp = sub_hilbert_spaces[to_spn];

 auto M = matrix<double>(to_sp.size(), from_sp.size());
 M() = 0;

 for (int i = 0; i<from_sp.size(); ++i) { // loop on all fock states of the blocks
  state<hilbert_space, double, true> from_s(full_hs);
  from_s(from_sp.get_fock_state(i)) = 1.0;

  auto to_s = op(from_s);
  auto proj_s = project<state<sub_hilbert_space, double, true>>(to_s,to_sp);

  // Count non-zero elements in this column
  foreach(proj_s, [&](int j, double ampl) { M(j,i) = ampl; nonzeros++; });
 }
 return {M,nonzeros};
}

// -------------------------------------------------------------------------------------------------

void sorted_spaces::complete_init1(many_body_op_t const& h_) {

 imperative_operator<hilbert_space, double, false> hamiltonian(h_, fops);

 /*
   Compute energy levels and eigenvectors of the local Hamiltonian
 */
 int n_subspaces = sub_hilbert_spaces.size();
 eigensystems.resize(n_subspaces);
 gs_energy = std::numeric_limits<double>::infinity();

 // Prepare the eigensystem in a temporary map to sort them by energy !
 std::map<std::pair<double, int>, eigensystem_t> eign_map;
 double energy_split = 1.e-10; // to split the eigenvalues which are numerically very close 
 for (int spn = 0; spn < n_subspaces; ++spn) {
  auto const& sp = sub_hilbert_spaces[spn];
  eigensystem_t eigensystem;
  if (_quantum_numbers.size()) eigensystem.quantum_numbers = _quantum_numbers[spn];

  state<sub_hilbert_space, double, false> i_state(sp);
  matrix<double> h_matrix(sp.size(), sp.size());

  for (int i = 0; i < sp.size(); ++i) {
   i_state.amplitudes()() = 0;
   i_state(i) = 1;
   auto f_state = hamiltonian(i_state);
   h_matrix(range(), i) = f_state.amplitudes();
  }

  auto eig = linalg::eigenelements(h_matrix);
  eigensystem.eigenvalues = eig.first;
  eigensystem.unitary_matrix = eig.second.transpose();
  gs_energy = std::min(gs_energy, eigensystem.eigenvalues[0]);

 /* eigensystem.eigenstates.reserve(sp.size());
  for (int e = 0; e < sp.size(); ++e) {
   eigensystem.eigenstates.emplace_back(sp);
   eigensystem.eigenstates.back().amplitudes() = h_matrix(e, range());
  }
*/
  eign_map.insert({{eigensystem.eigenvalues(0) + energy_split * spn, spn}, eigensystem});
 }

 // Reorder the block along their minimal energy
 {
  auto tmp = sub_hilbert_spaces;
  std::map<int, int> remap;
  int i = 0;
  for (auto const& x : eign_map) { // in order of min energy !
   eigensystems[i] = x.second;
   tmp[i] = sub_hilbert_spaces[x.first.second];
   tmp[i].set_index(i);
   remap[x.first.second] = i;
   ++i;
  }
  std::swap(tmp, sub_hilbert_spaces);
  auto remap_connection = [&](matrix<long>& connection) {
   auto c2 = connection;
   for (int n = 0; n < first_dim(connection); ++n)
    for (int i = 0; i < second_dim(connection); ++i) connection(n, remap[i]) = (c2(n, i) == -1 ? -1 : remap[c2(n, i)]);
  };
  remap_connection(creation_connection);
  remap_connection(annihilation_connection);
  // rematch the state which are NOT regular type !!
  //for (int spn = 0; spn < n_subspaces(); ++spn) {
  // for (auto& st : eigensystems[spn].eigenstates) st.set_hilbert(sub_hilbert_spaces[spn]);
  //}
 } // end reordering

 // Shift the ground state energy of the local Hamiltonian to zero.
 for (auto& eigensystem : eigensystems) eigensystem.eigenvalues() -= get_gs_energy();

 for (auto const& x : fops) {
  // get the operators and their index
  int n = x.linear_index;
  // n is guaranteed to be 0,1,2,3,  by the fundamental_operator_set class ... but it is very weird...
  // otherwise the push_back below is false.
  auto create = many_body_op_t::make_canonical(true, x.index);
  auto destroy = many_body_op_t::make_canonical(false, x.index);
  // construct their imperative counterpart
  imperative_operator<hilbert_space, double> op_c_dag(create, fops), op_c(destroy, fops);

  // Compute the matrices of c, c dagger in the diagonalization base of H_loc
  // first a lambda, since it is almost the same code for c and cdag
  auto make_c_mat = [&](int n, matrix<long> const& connection, imperative_operator<hilbert_space, double> c_op) {
   std::vector<matrix<double>> cmat(second_dim(connection));
   for (int B = 0; B < second_dim(connection); ++B) {
    auto Bp = connection(n, B);
    if (Bp == -1) continue;
    auto M = make_op_matrix(c_op, B, Bp);
    cmat[B] = eigensystems[Bp].unitary_matrix.transpose() * M.first * eigensystems[B].unitary_matrix;
   }
   return cmat;
  };

  // now execute code...
  c_matrices.push_back(make_c_mat(n, annihilation_connection, op_c));
  cdag_matrices.push_back(make_c_mat(n, creation_connection, op_c_dag));

 } // end of loop on operators

 _vacuum_index = -1;
 // get the position of the bare vacuum
 for (int bl = 0; bl < sub_hilbert_spaces.size(); ++bl) {
  if (sub_hilbert_spaces[bl].has_state(0)) {
   if (sub_hilbert_spaces[bl].size() != 1) TRIQS_RUNTIME_ERROR << "The bare vacuum is not in a block of size 1 !";
   _vacuum_index = bl;
   break;
  }
 }
 if (_vacuum_index < 0) TRIQS_RUNTIME_ERROR << "I did not find the bare vacuum !!";

 complete_init2();
 _quantum_numbers.clear();

}

// -----------------------------------------------------------------

void sorted_spaces::complete_init2() {
 _total_dim = 0;
 for (auto const& es : eigensystems) _total_dim += es.eigenvalues.size();

 // Calculate the index of the first eigenstate of each block
 first_eigstate_of_block.resize(_total_dim, 0);
 for (int bl = 1; bl < n_blocks(); ++bl) first_eigstate_of_block[bl] = first_eigstate_of_block[bl - 1] + get_block_dim(bl - 1);
}

// -----------------------------------------------------------------

std::string sorted_spaces::eigenstate_repr(int bl, int k) const {
 std::stringstream fs;
 auto const& U = eigensystems[bl].unitary_matrix;
 int n = first_dim(U);

 auto out = std::vector<std::string>{"0", "u", "d", "D"};

 for (int i = 0; i < n; ++i)
  if (std::abs(U(i, k)) > 1.e-10) {
   fs << "  ";
   std::bitset<16> v(sub_hilbert_spaces[bl].get_fock_state(i));
   for (int u = 3; u >= 0; --u) fs << out[v[2 * u] + 2 * v[2 * u + 1]];
   fs << "   " << U(i, k) << "\n";
  }
 return fs.str();
}


// -----------------------------------------------------------------
// Clean this vector vs std::vector
std::vector<std::tuple<vector<double>, std::vector<double>>> sorted_spaces::get_energy_and_quantum_numbers() const {
 decltype(get_energy_and_quantum_numbers()) R;
 for (auto const& es : eigensystems) R.push_back(std::make_tuple(es.eigenvalues, es.quantum_numbers));
 return R;
}

// -----------------------------------------------------------------

double sorted_spaces::partition_function(double beta) const {
 double z = 0;
 for (auto const& es : eigensystems)
  for (auto e : es.eigenvalues) z += std::exp(-beta * e);
 return z;
}

// -----------------------------------------------------------------

block_gf<imtime> sorted_spaces::atomic_gf(double beta, std::map<std::string, indices_t> const& gf_struct,
                                          int n_tau, std::vector<std::pair<int, int>> const& excluded_states) const {

 double z = partition_function(beta);

 std::vector<std::string> block_names;
 std::vector<gf<imtime>> gf_blocks;

 for (auto const& block : gf_struct) {
  block_names.push_back(block.first);
  int bl_size = block.second.size();
  auto g = gf<imtime>{{beta, Fermion, n_tau}, {bl_size,bl_size}};

  for(int inner_index1 = 0; inner_index1 < bl_size; ++inner_index1)
   for(int inner_index2 = 0; inner_index2 < bl_size; ++inner_index2){
    int n1 = fops[{block.first,block.second[inner_index1]}];    // linear_index of c
    int n2 = fops[{block.first,block.second[inner_index2]}];    // linear_index of c_dag

    for (int A = 0; A < n_blocks(); ++A) { // index of the A block. sum over all
     int B = creation_connection(n2,A);                  // index of the block connected to A by operator c_n
     if (B == -1) continue;                               // no matrix element
     if (annihilation_connection(n1,B) != A) continue;   //
     for (int ia = 0; ia < get_block_dim(A); ++ia)
      for (int ib = 0; ib < get_block_dim(B); ++ib){
       auto Ea = eigensystems[A].eigenvalues[ia];
       auto Eb = eigensystems[B].eigenvalues[ib];
       // exclusion
       auto be = excluded_states.begin(), ee = excluded_states.end();
       if ((std::find(be,ee, std::make_pair(A,ia)) != ee) || (std::find(be,ee, std::make_pair(B,ib)) != ee)) continue;
       for (auto tau : g.mesh())
        g[tau](inner_index1,inner_index2) += -cdag_matrices[n2][A](ib, ia) * c_matrices[n1][B](ia, ib) * std::exp(-(Eb - Ea) * tau - beta * Ea) / z;
      }
    }
  }
  g.singularity()(1) = 1.0;
  gf_blocks.push_back(g);
 }

 return make_block_gf(block_names, gf_blocks);
}

// -----------------------------------------------------------------

std::ostream& operator<<(std::ostream& os, sorted_spaces const& ss) {

 os << "Dimension of full Hilbert space: " << ss.dim() << std::endl;
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

block_matrix_t sorted_spaces::matrix_element(std::vector<int> P) {

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
  matrix<double> R(get_block_dim(bl), get_block_dim(bl));
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

block_matrix_t sorted_spaces::matrix_element(std::vector<std::pair<indices_t, indices_t>> const& P) {
 int L = fops.size();
 std::vector<int> Pi(L);
 for (int i = 0; i < L; ++i) Pi[i] = i;
 for (auto const& ii : P) Pi[fops[ii.first]] = fops[ii.second];
 return matrix_element(Pi);
}

// -----------------------------------------------------------------
std::pair<int, matrix<double>> sorted_spaces::matrix_element_of_monomial(many_body_op_t::monomial_t const& op_vec, int B) {

 matrix<double> m = triqs::arrays::make_unit_matrix<double>(get_block_dim(B));
 for (int i = op_vec.size() - 1; i >= 0; --i) {
  int ind = fops[op_vec[i].indices];
  m = m * (op_vec[i].dagger ? cdag_matrices[ind][B] : c_matrices[ind][B]);
  B = (op_vec[i].dagger ? creation_connection(ind, B) : annihilation_connection(ind, B));
  if (B == -1) break;
 }
 return {B, std::move(m)};
}

// -----------------------------------------------------------------

full_hilbert_space_state_t sorted_spaces::act(many_body_op_t const& op, full_hilbert_space_state_t const& st) {
 full_hilbert_space_state_t res(st.size());
 res() = 0;

 auto get_range = [this](int B) {
  int s = this->first_eigstate_of_block[B];
  return range{s, s + this->get_block_dim(B)};
 };

 for (auto const& x : op) {
  for (int bl = 0; bl < n_blocks(); ++bl) {
   auto b_m = matrix_element_of_monomial(x.monomial, bl);
   if (b_m.first != -1) res(get_range(b_m.first)) = x.coef * b_m.second * st(get_range(bl));
  }
 }
 return res;
}

// -----------------------------------------------------------------

namespace {
 // for two matrices
 double dot_product(matrix<double> const& a, matrix<double> const& b) {
  double r = 0;
  int dim1 = first_dim(a), dim2 = second_dim(a);
  if ((dim1 != second_dim(b)) || (dim2 != first_dim(b))) TRIQS_RUNTIME_ERROR << "dot_product of matrices : size mismatch";
  for (int i = 0; i < dim1; ++i)
   for (int j = 0; j < dim2; ++j) r += a(i, j) * b(j, i);
  return r;
 }
}

double sorted_spaces::average(block_matrix_t const& density_matrix, many_body_op_t const& op) {
 double r = 0;
 for (auto const& x : op)
  for (int bl = 0; bl < n_blocks(); ++bl) {
   auto b_m = matrix_element_of_monomial(x.monomial, bl);
   if (b_m.first != -1) r += x.coef * dot_product(b_m.second, density_matrix[bl]);
  }
 return r;
}

double sorted_spaces::average_on_projector(block_matrix_t const& density_matrix, full_hilbert_space_state_t const& psi) {
 double r = 0;
 for (int bl = 0; bl < n_blocks(); ++bl) {
  int d = get_block_dim(bl);
  for (int a = 0; a < d; ++a)
   for (int b = 0; b < d; ++b) r += psi(flatten_block_index(bl, a)) * density_matrix[bl](a, b) * psi(flatten_block_index(bl, b));
 }
 return r;
}
}
// move with the hpp file
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
std::string get_triqs_hdf5_data_scheme(sorted_spaces::eigensystem_t const&) { return "sorted_spaces::eigensystem_t"; }

// -----------------------------------------------------------------

void h5_write(h5::group fg, std::string const& name, sorted_spaces::eigensystem_t const& x) {
 auto gr = fg.create_group(name);
 h5_write(gr, "eigenvalues", x.eigenvalues);
 h5_write(gr, "unitary_matrix", x.unitary_matrix);
}

// -----------------------------------------------------------------
void h5_read(h5::group fg, std::string const& name, sorted_spaces::eigensystem_t& x) {
 auto gr = fg.open_group(name);
 h5_read(gr, "eigenvalues", x.eigenvalues);
 h5_read(gr, "unitary_matrix", x.unitary_matrix);
}

// -----------------------------------------------------------------
void h5_write(h5::group fg, std::string const& name, sorted_spaces const& x) {
 auto gr = fg.create_group(name);
 gr.write_triqs_hdf5_data_scheme(x);
 h5_write(gr, "creation_connection", x.creation_connection);
 h5_write(gr, "annihilation_connection", x.annihilation_connection);

 auto write_sparse = [&](std::string na, std::vector<std::vector<matrix<double>>> const& Mvv) {
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
 h5_write(gr, "vacuum_index", x._vacuum_index);
}

// -----------------------------------------------------------------
void h5_read(h5::group fg, std::string const& name, sorted_spaces& x) {
 auto gr = fg.open_group(name);
 h5_read(gr, "creation_connection", x.creation_connection);
 h5_read(gr, "annihilation_connection", x.annihilation_connection);

 auto read_sparse = [&](std::string na, std::vector<std::vector<matrix<double>>>& Mvv) {
  Mvv.resize(first_dim(x.creation_connection), std::vector<matrix<double>>(second_dim(x.creation_connection)));
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
 h5_read(gr, "vacuum_index", x._vacuum_index);
 x.complete_init2();
}
}
