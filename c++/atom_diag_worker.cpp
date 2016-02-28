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
#include "./atom_diag_worker.hpp"
#include <triqs/arrays/linalg/eigenelements.hpp>
#include <triqs/hilbert_space/space_partition.hpp>
#include <sstream>
#include <bitset>
#include <algorithm>
using namespace triqs::arrays;
using namespace triqs::hilbert_space;
using std::string;

namespace cthyb {

// Filter the fock states with a given number of particles
bool atom_diag_worker::fock_state_filter(fock_state_t s){
  auto c = std::bitset<64>(s).count();
  return ((c >= n_min) && (c <= n_max));
}

//------------------------------------------------------------------------------------

void atom_diag_worker::autopartition() { 

 fundamental_operator_set const& fops = hdiag->get_fops();
 many_body_op_t const& h = hdiag->get_h_atomic();
 
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
 hdiag->sub_hilbert_spaces.reserve(SP.n_subspaces());
 for (int n = 0; n < SP.n_subspaces(); ++n) hdiag->sub_hilbert_spaces.emplace_back(n);
  
 //FIXME foreach (SP, [&](fock_state_t s, int spn) { if (fock_state_filter(s)) hdiag->sub_hilbert_spaces[spn].add_fock_state(s); }) ;
 foreach (SP, [&](fock_state_t s, int spn) { hdiag->sub_hilbert_spaces[spn].add_fock_state(s); }) ;

 // Fill connections
 hdiag->creation_connection.resize(fops.size(), SP.n_subspaces());
 hdiag->annihilation_connection.resize(fops.size(), SP.n_subspaces());
 hdiag->creation_connection.as_array_view() = -1;
 hdiag->annihilation_connection.as_array_view() = -1;

 for (auto const& o : fops) {
  int n = o.linear_index;
  for (auto const& e : creation_melem[n]) {
   fock_state_t i, f;
   std::tie(i, f) = e.first;
   hdiag->creation_connection(n, SP.lookup_basis_state(i)) = SP.lookup_basis_state(f);
  }
  for (auto const& e : annihilation_melem[n]) {
   fock_state_t i, f;
   std::tie(i, f) = e.first;
   hdiag->annihilation_connection(n, SP.lookup_basis_state(i)) = SP.lookup_basis_state(f);
  }
 }
 complete();
}

// ************************************************************************************************

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

//-----------------------------

void atom_diag_worker::partition_with_qn(std::vector<many_body_op_t> const& qn_vector) {

 fundamental_operator_set const& fops = hdiag->get_fops();
 many_body_op_t const& h_ = hdiag->get_h_atomic();
 
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
   auto n_blocks = hdiag->sub_hilbert_spaces.size();
   hdiag->sub_hilbert_spaces.emplace_back(n_blocks); // a new sub_hilbert_space
   hdiag->quantum_numbers.push_back(qn);
   map_qn_n[qn] = n_blocks;
  }

  // add fock state to partial Hilbert space
  hdiag->sub_hilbert_spaces[map_qn_n[qn]].add_fock_state(fs);
 }

 // ----  now make the creation map -----

 auto creation_map = std::vector<std::vector<int>>(fops.size(), std::vector<int>(hdiag->sub_hilbert_spaces.size(), -1));
 auto annihilation_map = creation_map;

 // init the mapping tables
 hdiag->creation_connection.resize(fops.size(), hdiag->sub_hilbert_spaces.size());
 hdiag->annihilation_connection.resize(fops.size(), hdiag->sub_hilbert_spaces.size());
 hdiag->creation_connection.as_array_view() = -1;
 hdiag->annihilation_connection.as_array_view() = -1;

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
    auto origin = hdiag->sub_hilbert_spaces[map_qn_n[qn_before]].get_index();
    auto target = hdiag->sub_hilbert_spaces[map_qn_n[qn_after]].get_index();
    if (creation_map[n][origin] == -1)
     creation_map[n][origin] = target;
    else if (creation_map[n][origin] != target)
     TRIQS_RUNTIME_ERROR << "Internal Error, AtomicProblem, creation";
    hdiag->creation_connection(n,map_qn_n[qn_before]) = map_qn_n[qn_after];
   }

   // apply annihilation on state to figure quantum numbers
   qn_after = get_quantum_numbers(op_c(s));

   // insert in annihilation map checking whether it was already there
   if (!triqs::utility::is_zero(dot_product(op_c(s), op_c(s)))) {
    auto origin = hdiag->sub_hilbert_spaces[map_qn_n[qn_before]].get_index();
    auto target = hdiag->sub_hilbert_spaces[map_qn_n[qn_after]].get_index();
    if (annihilation_map[n][origin] == -1)
     annihilation_map[n][origin] = target;
    else if (annihilation_map[n][origin] != target)
     TRIQS_RUNTIME_ERROR << "Internal Error, AtomicProblem, annihilation";
    hdiag->annihilation_connection(n, map_qn_n[qn_before]) = map_qn_n[qn_after];
   }
  }
 }
 complete();
}

// -------------------------------------------------------------------------------------------------

matrix<double> atom_diag_worker::make_op_matrix(imperative_operator<hilbert_space> const& op, int from_spn,
                                                                int to_spn) const {

 fundamental_operator_set const& fops = hdiag->get_fops();
 hilbert_space full_hs(fops);
 auto const& from_sp = hdiag->sub_hilbert_spaces[from_spn];
 auto const& to_sp = hdiag->sub_hilbert_spaces[to_spn];

 auto M = matrix<double>(to_sp.size(), from_sp.size());
 M() = 0;

 for (int i = 0; i<from_sp.size(); ++i) { // loop on all fock states of the blocks
  state<hilbert_space, double, true> from_s(full_hs);
  from_s(from_sp.get_fock_state(i)) = 1.0;
  auto to_s = op(from_s);
  auto proj_s = project<state<sub_hilbert_space, double, true>>(to_s,to_sp);
  foreach(proj_s, [&](int j, double ampl) { M(j,i) = ampl;});
 }
 return M;
}

// ************************************************************************************************

void atom_diag_worker::complete() {

 fundamental_operator_set const& fops = hdiag->get_fops();
 many_body_op_t const& h_ = hdiag->get_h_atomic();

 imperative_operator<hilbert_space, double, false> hamiltonian(h_, fops);

 /*
   Compute energy levels and eigenvectors of the local Hamiltonian
 */
 int n_subspaces = hdiag->sub_hilbert_spaces.size();
 hdiag->eigensystems.resize(n_subspaces);
 hdiag->gs_energy = std::numeric_limits<double>::infinity();

 // Prepare the eigensystem in a temporary map to sort them by energy !
 std::map<std::pair<double, int>, atom_diag::eigensystem_t> eign_map;
 double energy_split = 1.e-10; // to split the eigenvalues which are numerically very close 
 for (int spn = 0; spn < n_subspaces; ++spn) {
  auto const& sp = hdiag->sub_hilbert_spaces[spn];
  atom_diag::eigensystem_t eigensystem;

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
  hdiag->gs_energy = std::min(hdiag->gs_energy, eigensystem.eigenvalues[0]);

//FIXME
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
  auto tmp = hdiag->sub_hilbert_spaces;
  std::map<int, int> remap;
  int i = 0;
  for (auto const& x : eign_map) { // in order of min energy !
   hdiag->eigensystems[i] = x.second;
   tmp[i] = hdiag->sub_hilbert_spaces[x.first.second];
   tmp[i].set_index(i);
   remap[x.first.second] = i;
   ++i;
  }
  std::swap(tmp, hdiag->sub_hilbert_spaces);
  auto remap_connection = [&](matrix<long>& connection) {
   auto c2 = connection;
   for (int n = 0; n < first_dim(connection); ++n)
    for (int i = 0; i < second_dim(connection); ++i) connection(n, remap[i]) = (c2(n, i) == -1 ? -1 : remap[c2(n, i)]);
  };
  remap_connection(hdiag->creation_connection);
  remap_connection(hdiag->annihilation_connection);
  // rematch the state which are NOT regular type !!
//FIXME
  //for (int spn = 0; spn < n_subspaces(); ++spn) {
  // for (auto& st : hdiag->eigensystems[spn].eigenstates) st.set_hilbert(hdiag->sub_hilbert_spaces[spn]);
  //}
 } // end reordering

 // Shift the ground state energy of the local Hamiltonian to zero.
 for (auto& eigensystem : hdiag->eigensystems) eigensystem.eigenvalues() -= hdiag->get_gs_energy();

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
    cmat[B] = hdiag->eigensystems[Bp].unitary_matrix.transpose() * M * hdiag->eigensystems[B].unitary_matrix;
   }
   return cmat;
  };

  // now execute code...
  hdiag->c_matrices.push_back(make_c_mat(n, hdiag->annihilation_connection, op_c));
  hdiag->cdag_matrices.push_back(make_c_mat(n, hdiag->creation_connection, op_c_dag));

 } // end of loop on operators

 hdiag->vacuum_block_index = -1;
 // get the position of the bare vacuum
 for (int bl = 0; bl < hdiag->sub_hilbert_spaces.size(); ++bl) {
  if (hdiag->sub_hilbert_spaces[bl].has_state(0)) {
   hdiag->vacuum_block_index = bl;
   hdiag->vacuum_inner_index = hdiag->sub_hilbert_spaces[bl].get_state_index(0);
   break;
  }
 }
 if (hdiag->vacuum_block_index < 0) TRIQS_RUNTIME_ERROR << "I did not find the bare vacuum !!";
}
}
