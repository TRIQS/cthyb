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
#include "sorted_spaces.hpp"
#include "./space_partition.hpp"
#include <fstream>
#include <triqs/arrays/linalg/eigenelements.hpp>

using namespace triqs::arrays;
using std::string;

namespace cthyb_matrix {

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

std::tuple<std::vector<sub_hilbert_space>, std::vector<std::vector<long>>, std::vector<std::vector<long>>>
autopartition(fundamental_operator_set const& fops, triqs::utility::many_body_operator<double> const& h) {

 imperative_operator<hilbert_space,false> hamiltonian(h,fops);
 hilbert_space full_hs(fops);
 state<hilbert_space,double,true> st(full_hs);

 using space_partition_t = space_partition<state<hilbert_space,double,true>,imperative_operator<hilbert_space,false>>;
 // Split the Hilbert space
 space_partition_t SP(st,hamiltonian,false);

 std::vector<space_partition_t::matrix_element_map_t> creation_melem(fops.n_operators());
 std::vector<space_partition_t::matrix_element_map_t> annihilation_melem(fops.n_operators());
 // Merge subspaces
 for (auto const& o : fops) {
  auto create = triqs::utility::many_body_operator<double>::make_canonical(true,o.index);
  auto destroy = triqs::utility::many_body_operator<double>::make_canonical(false,o.index);

  imperative_operator<hilbert_space> op_c_dag(create,fops), op_c(destroy,fops);

  int n = o.linear_index;
  std::tie(creation_melem[n],annihilation_melem[n]) = SP.merge_subspaces(op_c_dag,op_c,true);
 }

 // Fill subspaces
 std::vector<sub_hilbert_space> subspaces;
 subspaces.reserve(SP.n_subspaces());
 for(int n=0; n<SP.n_subspaces(); ++n) subspaces.emplace_back(n);

 foreach(SP,[&subspaces](fock_state_t s, int spn){ subspaces[spn].add_fock_state(s); });

 // Fill connections
 std::vector<std::vector<long>> creation_connection(fops.n_operators(),std::vector<long>(SP.n_subspaces(),-1));
 std::vector<std::vector<long>> annihilation_connection(fops.n_operators(),std::vector<long>(SP.n_subspaces(),-1));

 for (auto const& o : fops) {
  int n = o.linear_index;
  for(auto const& e : creation_melem[n]){
      fock_state_t i,f;
      std::tie(i,f) = e.first;
      creation_connection[n][SP.lookup_basis_state(i)] = SP.lookup_basis_state(f);
  }
  for(auto const& e : annihilation_melem[n]){
      fock_state_t i,f;
      std::tie(i,f) = e.first;
      annihilation_connection[n][SP.lookup_basis_state(i)] = SP.lookup_basis_state(f);
  }
 }

 return std::make_tuple(subspaces,creation_connection,annihilation_connection);
}

//-----------------------------

sorted_spaces::sorted_spaces(triqs::utility::many_body_operator<double> const& h_,
                             std::vector<triqs::utility::many_body_operator<double>> const& qn_vector,
                             fundamental_operator_set const& fops)
   : hamiltonian(h_, fops),
     creation_connection(fops.n_operators()),
     destruction_connection(fops.n_operators()),
     //creation_operators(fops.n_operators()),
     //destruction_operators(fops.n_operators()),
     fops(fops) {

 slice_hilbert_space_with_qn(h_, qn_vector, fops);
 complete_init(h_);
}

//-----------------------------

void sorted_spaces::slice_hilbert_space_with_qn(triqs::utility::many_body_operator<double> const& h_,
                                 std::vector<triqs::utility::many_body_operator<double>> const& qn_vector,
                                 fundamental_operator_set const& fops) {

 // hilbert spaces and quantum numbers
 std::map<std::vector<double>, int, lt_dbl> map_qn_n;
 using quantum_number_t = double;
 std::vector<std::vector<quantum_number_t>> quantum_numbers;

 // the full Hilbert space
 hilbert_space full_hs(fops);

 // The QN as operators : a vector of imperative operators for the quantum numbers
 std::vector<imperative_operator<hilbert_space>> qn_operators;
 for (auto& qn : qn_vector) qn_operators.emplace_back(qn, fops);

 // Helper function to get quantum numbers
 auto get_quantum_numbers = [&qn_operators](state<hilbert_space, double, true> const& s) {
  std::vector<quantum_number_t> qn;
  for (auto const& op : qn_operators) qn.push_back(dot_product(s, op(s)));
  return qn;
 };

 /*
   The first part consists in dividing the full Hilbert space
   into smaller subspaces using the quantum numbers
 */
 for (int r = 0; r < full_hs.dimension(); ++r) {

  // fock_state corresponding to r
  fock_state_t fs = full_hs.get_fock_state(r);

  // the state we'll act on
  state<hilbert_space, double, true> s(full_hs);
  s(r) = 1.0;

  // create the vector with the quantum numbers
  std::vector<quantum_number_t> qn = get_quantum_numbers(s);

  // if first time we meet these quantum numbers create partial Hilbert space
  if (map_qn_n.count(qn) == 0) {
   auto n_blocks = sub_hilbert_spaces.size();
   sub_hilbert_spaces.emplace_back(n_blocks); // a new sub_hilbert_space
   quantum_numbers.push_back(qn);
   map_qn_n[qn] = n_blocks;
  }

  // add fock state to partial Hilbert space
  sub_hilbert_spaces[map_qn_n[qn]].add_fock_state(fs);
 }


  /*
   In this second part we want to derive the partial Hilbert space
   mapping. Basically we want to know if we act on a partial Hilbert
   space with a creation (destruction) operator, in which other
   partial Hilbert space we end up.
 */

 auto creation_map = std::vector<std::vector<int>>(fops.n_operators(), std::vector<int>(sub_hilbert_spaces.size(), -1));
 auto destruction_map = creation_map;

 for (auto const& x : fops) {

  // get the operators and their index
  int n = x.linear_index;
  auto create = triqs::utility::many_body_operator<double>::make_canonical(true, x.index);
  auto destroy = triqs::utility::many_body_operator<double>::make_canonical(false, x.index);

  // construct their imperative counterpart
  imperative_operator<hilbert_space> op_c_dag(create, fops), op_c(destroy, fops);

  // to avoid declaring every time in the loop below
  std::vector<quantum_number_t> qn_before, qn_after;

  // these will be mapping tables
  creation_connection[n].resize(sub_hilbert_spaces.size(), -1);
  destruction_connection[n].resize(sub_hilbert_spaces.size(), -1);

  // now act on the state with the c, c_dag to see how quantum numbers change
  for (int r = 0; r < full_hs.dimension(); ++r) {

   // the state we'll act on and its quantum numbers
   state<hilbert_space, double, true> s(full_hs);
   s(r) = 1.0;
   qn_before = get_quantum_numbers(s);

   // apply creation on state to figure quantum numbers
   qn_after = get_quantum_numbers(op_c_dag(s));

   // insert in creation map checking whether it was already there
   if (dot_product(op_c_dag(s), op_c_dag(s)) > 1.e-10) {
    auto origin = sub_hilbert_spaces[map_qn_n[qn_before]].get_index();
    auto target = sub_hilbert_spaces[map_qn_n[qn_after]].get_index();
    if (creation_map[n][origin] == -1)
     creation_map[n][origin] = target;
    else if (creation_map[n][origin] != target)
     TRIQS_RUNTIME_ERROR << "Internal Error, Sorted Space, Creation";
    creation_connection[n][map_qn_n[qn_before]] = map_qn_n[qn_after];
   }

   // apply destruction on state to figure quantum numbers
   qn_after = get_quantum_numbers(op_c(s));

   // insert in destruction map checking whether it was already there
   if (dot_product(op_c(s), op_c(s)) > 1.e-10) {
    auto origin = sub_hilbert_spaces[map_qn_n[qn_before]].get_index();
    auto target = sub_hilbert_spaces[map_qn_n[qn_after]].get_index();
    if (destruction_map[n][origin] == -1)
     destruction_map[n][origin] = target;
    else if (destruction_map[n][origin] != target)
     TRIQS_RUNTIME_ERROR << "Internal Error, Sorted Space, Creation";
    destruction_connection[n][map_qn_n[qn_before]] = map_qn_n[qn_after];
   }
  }

  // IF WE CHOOSE TO KEEP, we MUST reorder ??
  // insert the creation and destruction operators in vectors. this is the fast version
  // of the operators because we explicitly use the map
  //creation_operators[n] = imperative_operator<sub_hilbert_space, true>(create, fops, creation_map[n], &sub_hilbert_spaces);
  //destruction_operators[n] = imperative_operator<sub_hilbert_space, true>(destroy, fops, destruction_map[n], &sub_hilbert_spaces);

 }
}

// -------------------------------------------------------------------------------------------------
 
 void sorted_spaces::complete_init(triqs::utility::many_body_operator<double> const& h_) { 

 /*
   Compute energy levels and eigenvectors of the local Hamiltonian
 */
 eigensystems.resize(n_subspaces());
 gs_energy = std::numeric_limits<double>::infinity();

 // Prepare the eigensystem in a temporary map to sort them by energy !
 std::map<std::pair<double, int>, eigensystem_t> eign_map;

 for (int spn = 0; spn < n_subspaces(); ++spn) {
  auto const& sp = subspace(spn);
  // auto& eigensystem = eigensystems[spn];
  eigensystem_t eigensystem;

  state<sub_hilbert_space, double, false> i_state(sp);
  matrix<double> h_matrix(sp.dimension(), sp.dimension());

  for (int i = 0; i < sp.dimension(); ++i) {
   i_state.amplitudes()() = 0;
   i_state(i) = 1;
   auto f_state = hamiltonian(i_state);
   h_matrix(range(), i) = f_state.amplitudes();
  }
  linalg::eigenelements_worker<matrix_view<double>, true> ew(h_matrix);

  ew.invoke();
  eigensystem.eigenvalues = ew.values();
  eigensystem.unitary_matrix = h_matrix.transpose();
  gs_energy = std::min(gs_energy, eigensystem.eigenvalues[0]);

  eigensystem.eigenstates.reserve(sp.dimension());
  for (int e = 0; e < sp.dimension(); ++e) {
   eigensystem.eigenstates.emplace_back(sp);
   eigensystem.eigenstates.back().amplitudes() = h_matrix(e, range());
  }

  eign_map.insert({{eigensystem.eigenvalues(0), spn}, eigensystem});
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
  //for (auto& x : map_qn_n) x.second = remap[x.second];
  auto remap_connection = [&](std::vector<std::vector<long>> & connection) {
   for (auto& cc : connection) {
    auto cc2 = cc;
    for (int i = 0; i < cc.size(); ++i) cc2[remap[i]] = remap[cc[i]];
    cc = cc2;
   }
  };
  remap_connection(creation_connection);
  remap_connection(destruction_connection);
  // rematch the state which are NOT regular type !!
  for (int spn = 0; spn < n_subspaces(); ++spn) {
   for (auto& st : eigensystems[spn].eigenstates) st.set_hilbert(sub_hilbert_spaces[spn]);
  }
 } // end reordering

 // Shift the ground state energy of the local Hamiltonian to zero.
 for (auto& eigensystem : eigensystems) eigensystem.eigenvalues() -= get_gs_energy();
 hamiltonian = imperative_operator<sub_hilbert_space, false>(h_ - get_gs_energy(), fops);

  // the full Hilbert space
 hilbert_space full_hs(fops);

 for (auto const& x : fops) {
  // get the operators and their index
  int n = x.linear_index;
  auto create = triqs::utility::many_body_operator<double>::make_canonical(true, x.index);
  auto destroy = triqs::utility::many_body_operator<double>::make_canonical(false, x.index);
  // construct their imperative counterpart
  imperative_operator<hilbert_space> op_c_dag(create, fops), op_c(destroy, fops);

   // Compute the matrices of c, c dagger in the diagonalization base of H_loc
  // first a lambda, since it is almost the same code for c and cdag
  auto make_c_mat = [&](std::vector<long> const& connection, imperative_operator<hilbert_space> c_op) {
   std::vector<triqs::arrays::matrix<double>> cmat(connection.size());
   for (int B = 0; B < connection.size(); ++B) {
    auto Bp = connection[B];
    if (Bp == -1) continue;
    auto M = matrix<double>(sub_hilbert_spaces[Bp].dimension(), sub_hilbert_spaces[B].dimension());
    M() = 0;
    // put the permutation matrix
    for (auto fs : sub_hilbert_spaces[B].get_all_fock_states()) { // loop on all fock states of the blocks
     state<hilbert_space, double, true> s(full_hs);
     s(fs) = 1.0;
     auto s2 = c_op(s);
     int nonzeros_in_s2 = 0;
     foreach(s2, [&](int i, double ampl){
         if (nonzeros_in_s2 >= 1) TRIQS_RUNTIME_ERROR << "Internal consistency error ";
         M(sub_hilbert_spaces[Bp].get_state_index(i/*full_hs.get_fock_state(i)==i*/), sub_hilbert_spaces[B].get_state_index(fs)) = ampl;
         nonzeros_in_s2++;
     });
    }
    cmat[B] = eigensystems[Bp].unitary_matrix.transpose() * M * eigensystems[B].unitary_matrix;
   }
   return cmat;
  };

  // now execute code...
  c_matrices.push_back(make_c_mat(destruction_connection[n], op_c));
  cdag_matrices.push_back(make_c_mat(creation_connection[n], op_c_dag));

 } // end of loop on operators

 std::ofstream("Diagonalization_atomic_pb") << *this;
}

// -----------------------------------------------------------------

double sorted_spaces::partition_function(double beta) const {
 // we should revert the order,
 double z = 0;
 for (auto const& es : eigensystems)
  for (auto e : es.eigenvalues) z += std::exp(-beta * e);
 return z;
}

// -----------------------------------------------------------------

block_gf<imtime> sorted_spaces::atomic_gf(double beta) const {

 double z = partition_function(beta);

 // std::cerr << "  Z = " << z << std::endl;
 auto n_times_pts = 100;
 std::vector<gf<imtime>> res;

 // this only valid if Gf is 1x1 matrices....
 for (auto const& x : fops) {

  int n = x.linear_index;
  auto g = gf<imtime>{{beta, Fermion, n_times_pts}, {1, 1}};

  for (int A = 0; A < sub_hilbert_spaces.size(); ++A) { // index of the A block. sum over all
   int B = creation_connection[n][A];                   // index of the block connected to A by operator c_n
   if (B == -1) continue;                               // no matrix element
   for (int ia = 0; ia < sub_hilbert_spaces[A].dimension(); ++ia)
    for (int ib = 0; ib < sub_hilbert_spaces[B].dimension(); ++ib) {
     auto Ea = eigensystems[A].eigenvalues[ia];
     auto Eb = eigensystems[B].eigenvalues[ib];
     for (auto tau : g.mesh())
      g[tau](0, 0) += -cdag_matrices[n][A](ib, ia) * c_matrices[n][B](ia, ib) * std::exp(-(Eb - Ea) * tau - beta * Ea) / z;
    }
  }
  res.push_back(g);
  // g /= - z; // nooverload ??
 }

 return make_block_gf(res);
}

// -----------------------------------------------------------------

std::ostream& operator<<(std::ostream& os, sorted_spaces const& ss) {

 os << "Number of blocks: " << ss.n_subspaces() << std::endl;
 for (int n = 0; n < ss.sub_hilbert_spaces.size(); ++n) {
  os << "Block " << n << ", ";
  os << "relative gs energy : " << ss.get_eigensystems()[n].eigenvalues[0] << std::endl;
  os << "size = " << ss.sub_hilbert_spaces[n].dimension() << std::endl;
  //os << "qn = ";
  //for (auto const& x : ss.quantum_numbers[n]) os << x << " ";
  os << std::endl;
  os << "index = " << ss.sub_hilbert_spaces[n].get_index() << std::endl;
  os << "-------------------------" << std::endl;
 }
 return os;
}
}
