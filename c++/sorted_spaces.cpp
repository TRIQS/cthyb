#include "sorted_spaces.hpp"

using namespace triqs::arrays;
using std::string;

namespace cthyb_krylov {

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

sorted_spaces::sorted_spaces(triqs::utility::many_body_operator<double> const& h_,
                             std::vector<triqs::utility::many_body_operator<double>> const& qn_vector,
                             fundamental_operator_set const& fops, std::vector<block_desc_t> const& block_structure)
   : hamiltonian(h_, fops),
     creation_operators(fops.n_operators()),
     destruction_operators(fops.n_operators()),
     creation_connection(fops.n_operators()),
     destruction_connection(fops.n_operators()) {

 std::map<indices_t, std::pair<int, int>> indices_to_ints;
 for (int bl = 0; bl < block_structure.size(); ++bl) {
  auto const& indices = block_structure[bl].indices;
  for (int i = 0; i < indices.size(); ++i) {
   indices_to_ints[indices[i]] = std::make_pair(bl, i);
  }
 }

 // hilbert spaces and quantum numbers
 std::map<std::vector<double>, int, lt_dbl> map_qn_n;

 // create the map int_pair_to_n : (int,int) --> int identifying operators
 for (auto x : fops) int_pair_to_n[indices_to_ints.at(x.index)] = x.linear_index;

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
   auto n_blocks = n_subspaces();
   sub_hilbert_spaces.emplace_back(n_blocks); // a new sub_hilbert_space
   quantum_numbers.push_back(qn);
   map_qn_n[qn] = n_blocks;
  }

  // add fock state to partial Hilbert space
  sub_hilbert_spaces[map_qn_n[qn]].add_fock_state(fs);
 }

 /* 
   Compute energy levels and eigenvectors of the local Hamiltonian
 */
 eigensystems.resize(n_subspaces());
 gs_energy = std::numeric_limits<double>::infinity();

 // Prepare the eigensystem in a temporary map to sort them by energy !
 std::map<std::pair<double,int>,eigensystem_t> eign_map;

 for (int spn = 0; spn < n_subspaces(); ++spn) {
  auto const& sp = subspace(spn);
  //auto& eigensystem = eigensystems[spn];
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
  for (auto& x : map_qn_n) {
   x.second = remap[x.second];
  }
  // rematch the state which are NOT regular type !!
  for (int spn = 0; spn < n_subspaces(); ++spn) {
   for (auto& st : eigensystems[spn].eigenstates) st.set_hilbert(sub_hilbert_spaces[spn]);
  }
 }

 // Shift the ground state energy of the local Hamiltonian to zero.
 for (auto& eigensystem : eigensystems) eigensystem.eigenvalues() -= get_gs_energy();
 hamiltonian = imperative_operator<sub_hilbert_space, false>(h_ - get_gs_energy(), fops);

 /*
   In this second part we want to derive the partial Hilbert space
   mapping. Basically we want to know if we act on a partial Hilbert
   space with a creation (destruction) operator, in which other
   partial Hilbert space we end up.
 */

 auto creation_map = std::vector<std::vector<int>>(fops.n_operators(), std::vector<int>(n_subspaces(), -1));
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
  creation_connection[n].resize(n_subspaces(), -1);
  destruction_connection[n].resize(n_subspaces(), -1);

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

  // insert the creation and destruction operators in vectors. this is the fast version
  // of the operators because we explicitly use the map
  creation_operators[n] = imperative_operator<sub_hilbert_space, true>(create, fops, creation_map[n], &sub_hilbert_spaces);
  destruction_operators[n] = imperative_operator<sub_hilbert_space, true>(destroy, fops, destruction_map[n], &sub_hilbert_spaces);
 }

}

// -----------------------------------------------------------------

std::ostream& operator<<(std::ostream& os, sorted_spaces const& ss) {

 os << "Number of blocks: " << ss.n_subspaces() << std::endl;
 for (int n = 0; n < ss.sub_hilbert_spaces.size(); ++n) {
  os << "Block " << n << ", ";
  os << "qn = ";
  for (auto const& x : ss.quantum_numbers[n]) os << x << " ";
  os << ", ";
  os << "index = " << ss.sub_hilbert_spaces[n].get_index() << std::endl ;
  os << "size = " << ss.sub_hilbert_spaces[n].dimension()<<std::endl ;
  os << " Relative gs energy : " << ss.get_eigensystems()[n].eigenvalues[0] << std::endl;
 }
 return os;
}

//----------------------------------------------------------------------
}
