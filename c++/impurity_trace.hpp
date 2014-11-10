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
#include "./configuration.hpp"
#include "./sorted_spaces.hpp"
#include "./solve_parameters.hpp"
#include "triqs/utility/rbt.hpp"
#include "triqs/statistics/histograms.hpp"

namespace cthyb {

/**
 */
class impurity_trace {
 public:

 int check_counter = 0; //DEBUG

 using trace_t = double;
 // using trace_t = std::complex<double>; TODO

 // construct from the config, the diagonalization of the loc Hamiltoninan, and parameters
 impurity_trace(configuration& c, sorted_spaces const& sosp, solve_parameters_t const& p);

 ~impurity_trace() { cancel_insert_impl(); } // in case of an exception, we need to remove any trial nodes before cleaning the tree!

 trace_t estimate(double p_yee = -1, double u_yee = 0);
 trace_t full_trace_over_estimator();

 private:
 // Various possible algorithms
 bool use_trace_estimator;
 enum class method_t {
  full_trace,
  estimate
 };
 method_t method;

 // ------- Configuration and h_loc data ----------------

 const configuration* config;            // config object does exist longer (temporally) than this object.
 const sorted_spaces* sosp;              // access to the diagonalization of h_loc
 int n_orbitals = sosp->n_c_operators(); //
 int n_blocks = sosp->n_subspaces();     //

 // ------------------ Cache data ----------------

 struct cache_t {
  double dtau_l = 0, dtau_r = 0;
  std::vector<int> block_table; // number of blocks limited to 2^15
  std::vector<arrays::matrix<double>> matrices;
  std::vector<double> matrix_lnorms; // -ln(norm(matrix))
  std::vector<bool> matrix_norm_valid;
  cache_t(int n_blocks) : block_table(n_blocks), matrix_lnorms(n_blocks), matrices(n_blocks), matrix_norm_valid(n_blocks) {}
 };

 struct node_data_t {
  op_desc op;
  cache_t cache;
  node_data_t(op_desc op, int n_blocks) : op(op), cache(n_blocks) {}
  void reset(op_desc op_new) { op = op_new; }
 };

 using rb_tree_t = rb_tree<time_pt, node_data_t, std::greater<time_pt>>;
 using node = rb_tree_t::node;

#ifdef EXT_DEBUG
 public:
#endif
 rb_tree_t tree;       // the red black tree and its nodes
 int n_modif;          // Analysis : number of nodes modified at the last change

 // ---------------- Cache machinery ----------------
 void update_cache();

 private:

 // The dimension of block b
 int get_block_dim(int b) const { return sosp->get_eigensystems()[b].eigenvalues.size(); }

 // the i-th eigenvalue of the block b
 double get_block_eigenval(int b, int i) const { return sosp->get_eigensystems()[b].eigenvalues[i]; }

 // the minimal eigenvalue of the block b
 double get_block_emin(int b) const { return get_block_eigenval(b, 0); }

 // node, block -> image of the block by n->op (the operator)
 int get_op_block_map(node n, int b) const {
  return sosp->fundamental_operator_connect(n->op.dagger, n->op.linear_index, b);
 }

 // the matrix of n->op, from block b to its image
 matrix<double> const& get_op_block_matrix(node n, int b) const {
  return sosp->fundamental_operator_matrix(n->op.dagger, n->op.linear_index, b);
 }

 // recursive function for tree traversal
 int compute_block_table(node n, int b);
 std::pair<int, double> compute_block_table_and_bound(node n, int b, double bound_threshold, bool use_threshold = true); //,double lnorm);
 std::pair<int, arrays::matrix<double>> compute_matrix(node n, int b);
 trace_t compute_trace(bool to_machine_precision, double p_yee, double u_yee);

 void update_cache_impl(node n);
 void update_dtau(node n);

 bool use_norm_of_matrices_in_cache = true; // When a matrix is computed in cache, its spectral radius replaces the norm estimate

 // integrity check
 void check_cache_integrity(bool print = false);
 void check_cache_integrity_one_node(node n, bool print);
 int check_one_block_table_linear(node n, int b, bool print); // compare block table to that of a linear method (ie. no tree)
 matrix<double> check_one_block_matrix_linear(node n, int b, bool print); // compare matrix to that of a linear method (ie. no tree)

 public:
 /*************************************************************************
  *  Ordinary binary search tree (BST) insertion of the trial nodes
  *************************************************************************/
 // We have a set of trial nodes, which we can glue, un-glue in the tree at will.
 // This avoids allocations.

 int tree_size = 0; // size of the tree +/- the added/deleted node

 // make a new detached black node
 std::shared_ptr<rb_tree_t::node_t> make_new_node() const {
  return std::make_shared<rb_tree_t::node_t>(time_pt{}, node_data_t{{}, n_blocks}, false, 1);
 }

 // a pool of trial nodes, ready to be glued in the tree. Max 4 (for possible double insertion)
 std::vector<std::shared_ptr<rb_tree_t::node_t>> trial_nodes = {make_new_node(), make_new_node(),
                                                                make_new_node(), make_new_node()};

 // for each inserted node, need to know {parent_of_node,child_is_left}
 std::vector<std::pair<node, bool>> inserted_nodes = {{nullptr, false}, {nullptr, false}, {nullptr, false}, {nullptr, false}};
 int trial_node_index = -1; // the index of the next available node in trial_nodes

 node try_insert_impl(node h, node n) { // implementation
  if (h == nullptr) return n;
  if (h->key == n->key) {
   // std::cerr << "insertion error "<< h->key <<  n->key;
   throw rbt_insert_error{};
  }
  auto smaller = tree.comparator()(n->key, h->key);
  if (smaller)
   h->left = try_insert_impl(h->left, n);
  else
   h->right = try_insert_impl(h->right, n);
  if (inserted_nodes[trial_node_index].first == nullptr) inserted_nodes[trial_node_index] = {h, smaller};
  h->modified = true;
  return h;
 }

 // unlink all glued trial nodes
 void cancel_insert_impl() {
  for (int i = 0; i <= trial_node_index; ++i) {
   auto& r = inserted_nodes[i];
   if (r.first != nullptr) (r.second ? r.first->left : r.first->right) = nullptr;
  }
  if (tree_size == trial_node_index + 1) tree.get_root() = nullptr;
 }

 /*************************************************************************
  * Node Insertion
  *************************************************************************/
 public:
 // Put a trial node at tau for operator op using an ordinary BST insertion (ie. not red black)
 void try_insert(time_pt const& tau, op_desc const& op) {
  if (trial_node_index > 3) TRIQS_RUNTIME_ERROR << "Error : more than 4 insertions ";
  auto& root = tree.get_root();
  inserted_nodes[++trial_node_index] = {nullptr, false};
  node n = trial_nodes[trial_node_index].get(); // get the next available node
  n->reset(tau, op);                            // change the time and op of the node
  root = try_insert_impl(root, n);              // insert it using a regular BST, no red black
  tree_size++;
 }

 // Remove all trial nodes from the tree
 void cancel_insert() {
  cancel_insert_impl();
  trial_node_index = -1;
  tree_size = tree.size();
  n_modif = tree.clear_modified();
  check_cache_integrity();
 }

 // confirm the insertion of the nodes, with red black balance
 void confirm_insert() {
  // remove BST inserted nodes
  cancel_insert_impl();
  // then reinsert the nodes used for real in rb tree
  for (int i = 0; i <= trial_node_index; ++i) {
   node n = trial_nodes[i].get();
   tree.insert(n->key, {n->op, n_blocks});
  }
  trial_node_index = -1;
  update_cache();
  tree_size = tree.size();
  n_modif = tree.clear_modified();
  check_cache_integrity();
 }

 /*************************************************************************
  * Node Removal
  *************************************************************************/
 private:
 std::vector<node> removed_ops;
 std::vector<time_pt> removed_key;

 public:
 // Find and mark as deleted the nth operator with fixed dagger and block_index
 // n=0 : first operator, n=1, second, etc...
 time_pt try_delete(int n, int block_index, bool dagger) noexcept {
  // traverse the tree, looking for the nth operator of the correct dagger, block_index
  int i = 0;
  node x = find_if(tree, [&](node no) {
   if (no->op.dagger == dagger && no->op.block_index == block_index) ++i;
   return i == n + 1;
  });
  removed_ops.push_back(x);      // store the node
  removed_key.push_back(x->key); // store the key
  tree.set_modified_from_root_to(x->key);
  x->delete_flag = true; // mark the node for deletion
  tree_size--;
  return x->key;
 }

 // Clean all the delete flags
 void cancel_delete() {
  for (auto& n : removed_ops) n->delete_flag = false;
  removed_ops.clear();
  removed_key.clear();
  tree_size = tree.size();
  n_modif = tree.clear_modified();
  check_cache_integrity();
 }

 // Confirm deletion: the nodes flagged for deletion are truly deleted
 void confirm_delete() {
  for (auto& k : removed_key) tree.delete_node(k); // CANNOT use the node here
  removed_ops.clear();
  removed_key.clear();
  update_cache();
  tree_size = tree.size();
  n_modif = tree.clear_modified();
  check_cache_integrity();
 }

 /*************************************************************************
  * Node shift (=insertion+deletion)
  *************************************************************************/

 // No try_shift implemented. Use combination of try_insert and try_delete instead.

 // Cancel the shift
 void cancel_shift() {

  // Inserted nodes
  cancel_insert_impl();
  trial_node_index = -1;

  // Deleted nodes
  for (auto& n : removed_ops) n->delete_flag = false;
  removed_ops.clear();
  removed_key.clear();

  tree_size = tree.size();
  n_modif = tree.clear_modified();
  check_cache_integrity();
 }

 // Confirm the shift of the node, with red black balance
 void confirm_shift() {

  // Inserted nodes
  //  first remove BST inserted nodes
  cancel_insert_impl();
  //  then reinsert the nodes used for real in rb tree
  for (int i = 0; i <= trial_node_index; ++i) {
   node n = trial_nodes[i].get();
   tree.insert(n->key, {n->op, n_blocks});
  }
  trial_node_index = -1;

  // Deleted nodes
  for (auto& k : removed_key) tree.delete_node(k); // CANNOT use the node here
  removed_ops.clear();
  removed_key.clear();

  // update cache only at the end
  update_cache();
  tree_size = tree.size();
  n_modif = tree.clear_modified();
  check_cache_integrity();
 }

 private:
 // ---------------- Histograms ----------------
 struct histograms_t {

  histograms_t(int n_subspaces) : n_subspaces(n_subspaces) {};
  int n_subspaces;

  // How many block non zero at root of the tree
  statistics::histogram n_block_at_root = {n_subspaces, "histo_n_block_at_root.dat"};

  // how many block kept after the truncation with the bound
  statistics::histogram n_block_kept = {n_subspaces, "histo_n_block_kept.dat"};

  // What is the dominant block in the trace computation ? Sorted by number or energy
  statistics::histogram dominant_block_bound = {n_subspaces, "histo_dominant_block_bound.dat"};
  statistics::histogram dominant_block_trace = {n_subspaces, "histo_dominant_block_trace.dat"};
  statistics::histogram_segment_bin dominant_block_energy_bound = {0, 100, 100, "histo_dominant_block_energy_bound.dat"};
  statistics::histogram_segment_bin dominant_block_energy_trace = {0, 100, 100, "histo_dominant_block_energy_trace.dat"};

  // Perturbation order, total and by color
  statistics::histogram opcount_total = {1000, "histo_opcount_total.dat"};
  std::vector<statistics::histogram> opcount;

  // Various ratios : trace/bound, trace/first term of the trace, etc..
  statistics::histogram_segment_bin trace_over_estimator = {0, 2, 100, "histo_trace_over_estimator.dat"};
  statistics::histogram_segment_bin trace_over_bound = {0, 1.5, 100, "histo_trace_over_bound.dat"};
  statistics::histogram_segment_bin trace_first_over_sec_term = {0, 1.0, 100, "histo_trace_first_over_sec_term.dat"};
  statistics::histogram_segment_bin trace_first_term_trace = {0, 1.0, 100, "histo_trace_first_term_trace.dat"};
 };
 std::unique_ptr<histograms_t> histo;
};
}
