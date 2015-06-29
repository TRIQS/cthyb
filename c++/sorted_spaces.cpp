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
#include "sorted_spaces.hpp"
#include <triqs/arrays/linalg/eigenelements.hpp>
#include <triqs/hilbert_space/space_partition.hpp>

using namespace triqs::arrays;
using namespace triqs::hilbert_space;
using std::string;

namespace cthyb {


//---------------------------------------------------------------------------------

void sorted_spaces::autopartition(fundamental_operator_set const& fops, many_body_op_t const& h) {

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

 foreach(SP, [this](fock_state_t s, int spn) { this->sub_hilbert_spaces[spn].add_fock_state(s); });

 // Fill connections
 creation_connection = std::vector<std::vector<long>>(fops.size(), std::vector<long>(SP.n_subspaces(), -1));
 annihilation_connection = std::vector<std::vector<long>>(fops.size(), std::vector<long>(SP.n_subspaces(), -1));

 for (auto const& o : fops) {
  int n = o.linear_index;
  for (auto const& e : creation_melem[n]) {
   fock_state_t i, f;
   std::tie(i, f) = e.first;
   creation_connection[n][SP.lookup_basis_state(i)] = SP.lookup_basis_state(f);
  }
  for (auto const& e : annihilation_melem[n]) {
   fock_state_t i, f;
   std::tie(i, f) = e.first;
   annihilation_connection[n][SP.lookup_basis_state(i)] = SP.lookup_basis_state(f);
  }
 }
}

//-----------------------------

sorted_spaces::sorted_spaces(many_body_op_t const& h_, std::vector<many_body_op_t> const& qn_vector,
                             fundamental_operator_set const& fops)
   : creation_connection(fops.size()), annihilation_connection(fops.size()), fops(fops), full_hs(fops) {
 slice_hilbert_space_with_qn(h_, qn_vector, fops);
 complete_init(h_);
}

//-----------------------------

sorted_spaces::sorted_spaces(many_body_op_t const& h_, fundamental_operator_set const& fops)
   : creation_connection(fops.size()), annihilation_connection(fops.size()), fops(fops), full_hs(fops) {
 autopartition(fops, h_);
 complete_init(h_);
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

 // hilbert spaces and quantum numbers
 std::map<std::vector<double>, int, lt_dbl> map_qn_n;
 using quantum_number_t = double;
 std::vector<std::vector<quantum_number_t>> quantum_numbers;

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
   quantum_numbers.push_back(qn);
   map_qn_n[qn] = n_blocks;
  }

  // add fock state to partial Hilbert space
  sub_hilbert_spaces[map_qn_n[qn]].add_fock_state(fs);
 }

 // ----  now make the creation map -----

 auto creation_map = std::vector<std::vector<int>>(fops.size(), std::vector<int>(sub_hilbert_spaces.size(), -1));
 auto annihilation_map = creation_map;

 for (auto const& x : fops) {

  // get the operators and their index
  int n = x.linear_index;
  auto create = many_body_op_t::make_canonical(true, x.index);
  auto destroy = many_body_op_t::make_canonical(false, x.index);

  // construct their imperative counterpart
  imperative_operator<hilbert_space, double> op_c_dag(create, fops), op_c(destroy, fops);

  // to avoid declaring every time in the loop below
  std::vector<quantum_number_t> qn_before, qn_after;

  // these will be mapping tables
  creation_connection[n].resize(sub_hilbert_spaces.size(), -1);
  annihilation_connection[n].resize(sub_hilbert_spaces.size(), -1);

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
    creation_connection[n][map_qn_n[qn_before]] = map_qn_n[qn_after];
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
    annihilation_connection[n][map_qn_n[qn_before]] = map_qn_n[qn_after];
   }
  }
 }
}

// -------------------------------------------------------------------------------------------------
std::pair<matrix<double>,int> sorted_spaces::make_op_matrix(imperative_operator<hilbert_space> const& op, int from_spn, int to_spn) const {

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

void sorted_spaces::complete_init(many_body_op_t const& h_) {

 imperative_operator<hilbert_space, double, false> hamiltonian(h_, fops);

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

  eigensystem.eigenstates.reserve(sp.size());
  for (int e = 0; e < sp.size(); ++e) {
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
  auto remap_connection = [&](std::vector<std::vector<long>>& connection) {
   for (auto& cc : connection) {
    auto cc2 = cc;
    for (int i = 0; i < cc.size(); ++i) cc2[remap[i]] = (cc[i] == -1 ? -1 : remap[cc[i]]);
    cc = cc2;
   }
  };
  remap_connection(creation_connection);
  remap_connection(annihilation_connection);
  // rematch the state which are NOT regular type !!
  for (int spn = 0; spn < n_subspaces(); ++spn) {
   for (auto& st : eigensystems[spn].eigenstates) st.set_hilbert(sub_hilbert_spaces[spn]);
  }
 } // end reordering

 // Shift the ground state energy of the local Hamiltonian to zero.
 for (auto& eigensystem : eigensystems) eigensystem.eigenvalues() -= get_gs_energy();
 // hamiltonian = imperative_operator<sub_hilbert_space, false>(h_ - get_gs_energy(), fops);

 for (auto const& x : fops) {
  // get the operators and their index
  int n = x.linear_index;
  auto create = many_body_op_t::make_canonical(true, x.index);
  auto destroy = many_body_op_t::make_canonical(false, x.index);
  // construct their imperative counterpart
  imperative_operator<hilbert_space, double> op_c_dag(create, fops), op_c(destroy, fops);

  // Compute the matrices of c, c dagger in the diagonalization base of H_loc
  // first a lambda, since it is almost the same code for c and cdag
  auto make_c_mat = [&](std::vector<long> const& connection, imperative_operator<hilbert_space, double> c_op) {
   std::vector<matrix<double>> cmat(connection.size());
   for (int B = 0; B < connection.size(); ++B) {
    auto Bp = connection[B];
    if (Bp == -1) continue;
    auto M = make_op_matrix(c_op,B,Bp);
    cmat[B] = eigensystems[Bp].unitary_matrix.transpose() * M.first * eigensystems[B].unitary_matrix;
   }
   return cmat;
  };

  // now execute code...
  c_matrices.push_back(make_c_mat(annihilation_connection[n], op_c));
  cdag_matrices.push_back(make_c_mat(creation_connection[n], op_c_dag));

 } // end of loop on operators
}

// -----------------------------------------------------------------

double sorted_spaces::partition_function(double beta) const {
 double z = 0;
 for (auto const& es : eigensystems)
  for (auto e : es.eigenvalues) z += std::exp(-beta * e);
 return z;
}

// -----------------------------------------------------------------

block_gf<imtime> sorted_spaces::atomic_gf(double beta, std::map<std::string,many_body_op_t::indices_t> const& gf_struct, int n_tau) const {


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

    for (int A = 0; A < sub_hilbert_spaces.size(); ++A) { // index of the A block. sum over all
     int B = creation_connection[n2][A];                  // index of the block connected to A by operator c_n
     if (B == -1) continue;                               // no matrix element
     if (annihilation_connection[n1][B] != A) continue;   //
     for (int ia = 0; ia < sub_hilbert_spaces[A].size(); ++ia)
      for (int ib = 0; ib < sub_hilbert_spaces[B].size(); ++ib){
       auto Ea = eigensystems[A].eigenvalues[ia];
       auto Eb = eigensystems[B].eigenvalues[ib];
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

std::map<std::string,std::vector<double>> sorted_spaces::atomic_observables(std::map<std::string,many_body_op_t> const& obs_map) const {

 std::map<std::string,std::vector<double>> results;
 std::map<std::string,imperative_operator<hilbert_space, double>> imp_ops;
 for(auto const& op : obs_map) {
  imp_ops.insert({op.first,{op.second,fops}});
  auto & res = results[op.first];
  res.reserve(full_hs.size());
 }

 for(auto const& eigs: eigensystems) {
  for(auto const& estate : eigs.eigenstates) {
   // eigenstate projected to the full Hilbert space
   auto proj_estate = project<state<hilbert_space, double, true>>(estate,full_hs);
   for(auto const& op : obs_map) results[op.first].push_back(dot_product(proj_estate,imp_ops[op.first](proj_estate)));
  }
 }

 return results;
}

// -----------------------------------------------------------------

std::ostream& operator<<(std::ostream& os, sorted_spaces const& ss) {

 os << "Dimension of full Hilbert space: " << ss.space().size() << std::endl;
 os << "Number of blocks: " << ss.n_subspaces() << std::endl;
 for (int n = 0; n < ss.sub_hilbert_spaces.size(); ++n) {
  os << "Block " << n << ", ";
  os << "relative gs energy : " << ss.get_eigensystems()[n].eigenvalues[0] << std::endl;
  os << "size = " << ss.sub_hilbert_spaces[n].size() << std::endl;
  // os << "qn = ";
  // for (auto const& x : ss.quantum_numbers[n]) os << x << " ";
  os << std::endl;
  os << "index = " << ss.sub_hilbert_spaces[n].get_index() << std::endl;
  // os << "             --------------" << std::endl;
  // for (auto const & cc : ss.creation_connection)  for (auto const & x : cc)  os << x << std::endl;
  // for (auto const& cc : ss.annihilation_connection) for (auto const & x : cc)  os << x << std::endl;
  os << "-------------------------" << std::endl;
 }
 return os;
}
}
