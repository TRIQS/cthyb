#pragma once
#include "configuration.hpp"
#include "sorted_spaces.hpp"
#include "rbt.hpp"
#include <triqs/parameters.hpp>
#include "triqs/statistics/histograms.hpp"

namespace cthyb_krylov {

/**
 * A worker that computes the trace using krylov method, for a given configuration.
 * Has to live longer than the configuration...
 */
class atomic_correlators_worker {
 public:
 using trace_t = double;
 // using trace_t = std::complex<double>;

 // construct from the config, the diagonalization of the loc Hamiltoninan, and parameters
 atomic_correlators_worker(configuration& c, sorted_spaces const& sosp, utility::parameters const& p);

 ~atomic_correlators_worker() { unlink_trial_nodes(); }
 // in case of an exception, we need to remove the trial node before cleaning the tree !!

 trace_t estimate();
 trace_t full_trace_over_estimator();

 private:
 // Various possible algorithms
 enum class method_t {
  FullTrace,
  EstimateWithBounds,
  EstimateTruncEps
 };
 method_t method;

 //------- data ----------------

 const configuration* config;            // must exists longer than this object.
 const sorted_spaces* sosp;              // access to the diagonalization of Hloc
 int n_orbitals = sosp->n_c_operators(); // 
 int n_blocks = sosp->n_subspaces();     //

 // ------------------ Cache data ----------------

 struct cache_t {
  double dtl = 0, dtr = 0;
  std::vector<int> block_table; //number of blocks limited to 2^15
  std::vector<double> matrix_lnorms; // - ln (norm(matrix)) 
  std::vector<arrays::matrix<double>> matrices;
  std::vector<bool> matrix_norm_are_valid;
  cache_t(int n_blocks) : block_table(n_blocks), matrix_lnorms(n_blocks), matrices(n_blocks), matrix_norm_are_valid(n_blocks) {}
 };

 struct node_data_t {
  op_desc op;
  cache_t cache;
  node_data_t(op_desc op, int n_blocks) : op(op), cache(n_blocks){}
  void reset(op_desc op1) { op = op1; }
 };

 using rb_tree_t = rb_tree<time_pt, node_data_t, std::greater<time_pt>>;
 using node = rb_tree_t::node;

#ifdef EXT_DEBUG 
 public:
#endif
 rb_tree_t tree; // the red black tree and its nodes
 int n_modif; // number of node modified at the last change
 void cache_update();

 /*************************************************************************
  *  Ordinary BST insertion of the trial nodes
  *************************************************************************/
 // we have a set of trial nodes, which we can glue, un-glue in the tree at will
 // avoids allocations.

 int tree_size = 0; // size of the tree +/- the added/deleted node

 // make a new detached black node
 std::shared_ptr<rb_tree_t::node_t> make_new_node() const {
  return std::make_shared<rb_tree_t::node_t>(time_pt{}, node_data_t{{}, n_blocks}, false, 1);
 }

 // a pool of trial nodes, ready to be glued in the tree. Max 4 (for possible double insertion)
 std::vector<std::shared_ptr<rb_tree_t::node_t>> trial_nodes = {make_new_node(), make_new_node(),
                                                                make_new_node(), make_new_node()};
 // {parent_of_node,child_is_left}
 std::vector<std::pair<node, bool>> inserted_nodes = {{nullptr, false}, {nullptr, false}, {nullptr, false}, {nullptr, false}};
 int trial_node_index = -1; // the index of the next available node in trial_nodes

 node bst_ordinary_insert(node h, node n) { // implementation
  if (h == nullptr) return n;
  if (h->key == n->key) throw rbt_insert_error{};
  auto smaller = tree.comparator()(n->key, h->key);
  if (smaller)
   h->left = bst_ordinary_insert(h->left, n);
  else
   h->right = bst_ordinary_insert(h->right, n);
  if (inserted_nodes[trial_node_index].first == nullptr) inserted_nodes[trial_node_index] = {h, smaller};
  h->modified = true;
  return h;
 }

 // unlink all glued trial nodes
 void unlink_trial_nodes() {
  for (int i = 0; i <= trial_node_index; ++i) {
   auto & r = inserted_nodes[i];
   if (r.first != nullptr) (r.second ? r.first->left : r.first->right) = nullptr;
  }
  if (tree_size == trial_node_index + 1) tree.get_root() = nullptr;
 }

 public:
 // Put a trial node at tau, op, using an ordinary BST insertion (no red black)
 void trial_node_insert(time_pt const& tau, op_desc const& op) {
  if (trial_node_index > 3) TRIQS_RUNTIME_ERROR << "Error : more than 4 insertions ";
  auto & root = tree.get_root();
  inserted_nodes[++trial_node_index] = {nullptr, false};
  node n = trial_nodes[trial_node_index].get(); // get the next available node
  n->reset(tau, op);                              // change the time and op of the node
  root = bst_ordinary_insert(root, n);            // insert it using a regular BST, no red black
  tree_size++;
 
  //recompute the dt
  // n->cache.dtr = (n->right ? double(n->key - tree.min_key(n->right)) : 0);
  // n->cache.dtl = (n->left ? double(tree.max_key(n->left) - n->key) : 0);
  }

 // Remove all trial nodes from the tree
 void trial_node_uninsert() {
  unlink_trial_nodes();
  trial_node_index = -1;
  tree_size = tree.size();
  n_modif = tree.clear_modified();
  check_cache_integrity(false);
 }

 // confirm the insertion of the nodes, with red black balance
 void confirm_trial_node_insertion() {
  unlink_trial_nodes();
  // then reinsert the nodes used for real in rb tree
  for (int i = 0; i <= trial_node_index; ++i) {
   node n = trial_nodes[i].get();
   tree.insert(n->key, {n->op, n_blocks});
  }
  trial_node_index = -1;
  cache_update();
  tree_size = tree.size();
  n_modif = tree.clear_modified();
  check_cache_integrity(false);
 }

 /*************************************************************************
  *  soft delete
  *************************************************************************/
 private:
 std::vector<node> removed_ops;

 public:
 // Find and soft_delete the n-th operator with fixed dagger and block_index
 // n=0 : first operator, n=1, second, etc...
 time_pt soft_delete_n_th_operator(int n, int block_index, bool dagger) {
  // traverse the tree, looking for the nth operator of the correct dagger, block_index
  int i = 0;
  node x = find_if(tree, [&](node no) {
   if (no->op.dagger == dagger && no->op.block_index == block_index) ++i;
   return i == n + 1;
  });
  removed_ops.push_back(x); // store the node
  tree.set_modified_from_root_to(x->key);
  x->soft_deleted = true;   // mark the node for deletion
  tree_size--;
  return x->key;
 }

 /// Clean all the soft deleted flags
 void clean_soft_delete() {
 for (auto& n : removed_ops) n->soft_deleted = false;
  removed_ops.clear();
  tree_size = tree.size();
  n_modif = tree.clear_modified();
  check_cache_integrity(false);
 }

 /// Confirm deletion : the soft deleted flagged node are truly deleted
 void confirm_soft_delete() {
  for (auto& n : removed_ops) tree.delete_node(n->key);
  removed_ops.clear();
  cache_update();
  tree_size = tree.size();
  n_modif = tree.clear_modified();
  check_cache_integrity(false);
 }

 private:
 // ---------------- Cache machinery ----------------

 // node, block -> image of the block by n->op (the operator)
 int get_op_block_map(node n, int b) const {
  return sosp->fundamental_operator_connect_from_linear_index(n->op.dagger, n->op.linear_index, b);
 }

 // The dimension of block b
 int get_block_dim(int b) const { return sosp->get_eigensystems()[b].eigenvalues.size(); }
 
 // the i-th eigenvalue of the block b
 double get_block_eigenval(int b, int i) const { return sosp->get_eigensystems()[b].eigenvalues[i]; }

 // the minimal eigenvalue of the block i
 double get_block_emin(int b) const { return get_block_eigenval(b, 0); }

 // the matrix of n->op, from block b to its image
 matrix<double> const& get_op_block_matrix(node n, int b) const {
  return sosp->fundamental_operator_matrix_from_linear_index(n->op.dagger, n->op.linear_index, b);
 }

 //recursive function for tree traversal
 int compute_block_table(node n, int b);
 std::pair<int, double> compute_block_table_and_bound(node n, int b, double bound_threshold);//,double lnorm);
 std::pair<int, arrays::matrix<double>> compute_matrix(node n, int b);

 // integrity check
 void check_cache_integrity(bool temporary, bool print=false);
 void check_cache_integrity_one_node(node n, bool temporary, bool print);
 int check_one_block_table_linear(node n, int b, bool print);
 matrix<double> check_one_block_matrix_linear(node n, int b, bool print);

 trace_t compute_trace(double epsilon = 1.e-15, bool estimator_only = false);
 trace_t last_estimate;
 
 //void update_cache();
 void update_cache_impl (node n);
 void update_dt(node n);

 double trace_epsilon_estimator() { return 1.e-1; }

 int block_dominant = 0;

 // ---------------- Histograms ----------------
 bool make_histograms;                       // Do we make the Histograms ?
 bool use_truncation;                        //
 //bool use_norm_of_matrices_in_cache = true; // When a matrix is computed in cache, its spectral radius replaces the norm estimate
 bool use_norm_of_matrices_in_cache = false; // When a matrix is computed in cache, its spectral radius replaces the norm estimate
 //bool use_only_first_term_in_trace = true; // Use only the first term in the trace 
 bool use_only_first_term_in_trace = false; //
 double double_max = std::numeric_limits<double>::max(); // easier to read

 //data for histo
 int n_call=0;

 // histograms
 statistics::histogram histo_dims_1x1 = {2, "histo_dims_1x1.dat"};
 statistics::histogram histo_nblock_at_root = {sosp->n_subspaces(), "histo_nblock_at_root.dat"};
 statistics::histogram histo_n_block_kept = {sosp->n_subspaces(), "histo_n_block_kept.dat"};
 statistics::histogram histo_n_block_examined_first_pass = {sosp->n_subspaces(), "histo_n_block_examined_first_pass.dat"};
 statistics::histogram histo_dominant_block_bound = {sosp->n_subspaces(), "histo_dominant_block_bound.dat"};
 statistics::histogram histo_dominant_block_trace = {sosp->n_subspaces(), "histo_dominant_block_trace.dat"};
 statistics::histogram histo_opcount_total = {1000, "histo_opcount_total.dat"};
 statistics::histogram histo_appel = {sosp->n_subspaces(), "histo_appel.dat"};
 statistics::histogram histo_appel2 = {50, "histo_appel2.dat"};
 statistics::histogram histo_appel3 = {50, "histo_appel3.dat"};
 
 statistics::histogram histo_cut1 = {sosp->n_subspaces(), "histo_cut1.dat"};
 statistics::histogram histo_cut2 = {sosp->n_subspaces(), "histo_cut2.dat"};

 statistics::histogram_segment_bin histo_trace_over_bound = {0, 1.5, 100, "hist_trace_over_bound.dat"};
 //statistics::histogram_segment_bin histo_perte_prod = {0, 1.5, 100, "histo_perte_prod.dat"};
 statistics::histogram_segment_bin histo_trace_first_over_sec_term = {0, 1.0, 100, "histo_trace_first_over_sec_term.dat"};
 statistics::histogram_segment_bin histo_trace_first_term_trace = {0, 1.0, 100, "histo_trace_first_term_trace.dat"};
 statistics::histogram_segment_bin histo_dominant_block_energy_bound = {0, 100, 100, "histo_dominant_block_energy_bound.dat"};
 statistics::histogram_segment_bin histo_dominant_block_energy_trace = {0, 100, 100, "histo_dominant_block_energy_trace.dat"};
 statistics::histogram_segment_bin histo_cut2_energy = {0,100,100, "histo_cut2_energy.dat"};
 
 std::vector<statistics::histogram> histo_opcount;
 //std::vector<statistics::histogram> histo_n_blocks_after_steps;
};
}
