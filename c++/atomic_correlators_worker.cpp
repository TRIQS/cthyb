#include "atomic_correlators_worker.hpp"
#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>
#include <algorithm>
#include <limits>

//#define CHECK_ALL

#ifdef CHECK_ALL
#define CHECK_CACHE
#define CHECK_AGAINST_LINEAR_COMPUTATION
#endif

// --------------- Computation of the matrix norm --> move into triqs::arrays ------------------------

double norm_induced_2(matrix<double> const& A) {

 // WORKAROUND BUG !!
 matrix<double> AA = A.transpose();
 auto M = AA * A;
 // auto M = A.transpose()* A;
 triqs::arrays::linalg::eigenelements_worker<matrix_view<double>, true> w(M());
 w.invoke();
 auto const& Es = w.values();
 return std::sqrt(Es(first_dim(Es) - 1)); // ordered is guaranteed by lapack
}


namespace cthyb_krylov {

atomic_correlators_worker::atomic_correlators_worker(configuration& c, sorted_spaces const& sosp_, utility::parameters const& p)
   : config(&c), sosp(&sosp_) {

 // Taking parameters from the inputs
 use_truncation = p["use_truncation"];

 std::string ms = p["trace_estimator"];
 try {
  method = std::map<std::string, method_t> {
   { "FullTrace", method_t::FullTrace }
   , {"EstimateWithBounds", method_t::EstimateWithBounds}, { "EstimateTruncEps", method_t::EstimateTruncEps }
  }
  .at(ms);
 }
 catch (...) {
  TRIQS_RUNTIME_ERROR << "Trace method " << ms << " not recognized.";
 }

 make_histograms = p["make_path_histograms"];
 if (make_histograms) {
 }
}
//------------------------------------------------------------------------------

atomic_correlators_worker::trace_t atomic_correlators_worker::estimate() {
 if (method == method_t::FullTrace) last_estimate = compute_trace(1.e-15, false);
 if (method == method_t::EstimateWithBounds) last_estimate = compute_trace(1.e-15, true);
 if (method == method_t::EstimateTruncEps) last_estimate = compute_trace(trace_epsilon_estimator(), false);
 return last_estimate;
}

//------------------------------------------------------------------------------

atomic_correlators_worker::trace_t atomic_correlators_worker::full_trace_over_estimator() {
 trace_t r = 1;
 if (method == method_t::EstimateWithBounds) r = compute_trace(1.e-15, false) / last_estimate;
 if (method == method_t::EstimateTruncEps) r = compute_trace(1.e-15, false) / last_estimate;
 if (make_histograms) histos["FullTrace_over_Estimate"] << std::abs(r);
 return r;
}

// ------------------- computation block table only -------------
/*
// for subtree at node n, return (B', bound)
int atomic_correlators_worker::compute_block_table(node n, int b) {
 if (b == -1) return -1;
 if (n == nullptr) return b; // identity
 if (!n->modified) return n->cache.block_table[b];
 int b1 = compute_block_table(n->right, b);
 if (b1 == -1) return -1;
 int b2 = (n->soft_deleted ? b1 : get_op_block_map(n, b1));
 if (b2 == -1) return -1;
 int b3 = compute_block_table(n->left, b2);
 return b3;
}
*/
// ------------------- computation block table and bounds -------------

// for subtree at node n, return (B', bound)
std::pair<int, double> atomic_correlators_worker::compute_block_table_and_bound(node n, int b, double lnorm_threshold) {

 if (b == -1) return {-1, 0};
 if (n == nullptr) return {b, 0}; // identity
 if (!n->modified) return {n->cache.block_table[b], n->cache.matrix_lnorms[b]};
 
 double lnorm1 = 0, lnorm3 = 0;
 double a_r = 0, a_l = 0;

 int b1 = b;
 if (n->right) {
  std::tie(b1, lnorm1) = compute_block_table_and_bound(n->right, b, lnorm_threshold);
  if (b1 == -1) return {-1, 0};
  a_r = double(n->key - tree.min_key(n->right)) * get_block_emin(b1);
 }

 int b2 = (n->soft_deleted ? b1 : get_op_block_map(n, b1));
 if (b2 == -1) return {-1, 0};

 int b3 = b2;
 if (n->left) {
  std::tie(b3, lnorm3) = compute_block_table_and_bound(n->left, b2, lnorm_threshold);
  if (b3 == -1) return {-1, 0};
  a_l = double(tree.max_key(n->left) - n->key) * get_block_emin(b2);
 }

 auto lnorm = lnorm3 + a_l + a_r + lnorm1;
 if ((lnorm_threshold > 0) && (lnorm > lnorm_threshold)) return {-1, 0};
 return {b3, lnorm};
}

// --------------- Computation of the matrix ------------------------------

std::pair<int, arrays::matrix<double>> atomic_correlators_worker::compute_matrix(node n, int b) {

 if (b == -1) return {-1, {}};
 if (n == nullptr) return {b, {}};
 if (!n->modified && n->cache.matrix_norm_are_valid[b]) return {n->cache.block_table[b], n->cache.matrices[b]};
 bool updating = (!n->modified && !n->cache.matrix_norm_are_valid[b]);

 double dtau_l = 0, dtau_r = 0;
 auto _ = arrays::range();

 auto r = compute_matrix(n->right, b);
 int b1 = r.first;
 if (b1 == -1) return {-1, {}};

 int b2 = (n->soft_deleted ? b1 : get_op_block_map(n, b1));
 if (b2 == -1) return {-1, {}};

 matrix<double> M = (!n->soft_deleted ? get_op_block_matrix(n, b1) : make_unit_matrix<double>(get_block_dim(b1)));

 if (n->right) { // M <- M* exp * r[b]
  dtau_r = double(n->key - tree.min_key(n->right));
  auto dim = second_dim(M); // same as get_block_dim(b2);
  for (int i = 0; i < dim; ++i) M(_, i) *= std::exp(-dtau_r * get_block_eigenval(b1, i));
  M = M * r.second; // optimise lapack call ?
 }

 int b3 = b2;
 if (n->left) { // M <- l[b] * exp * M
  auto l = compute_matrix(n->left, b2);
  b3 = l.first;
  if (b3 == -1) return {-1, {}};
  dtau_l = double(tree.max_key(n->left) - n->key);
  auto dim = first_dim(M); // same as get_block_dim(b1);
  for (int i = 0; i < dim; ++i) M(i, _) *= std::exp(-dtau_l * get_block_eigenval(b2, i));
  M = l.second * M;
 }

 if (updating) {
  n->cache.matrices[b] = M;
  n->cache.matrix_norm_are_valid[b] = true;

  // improve the norm
  if (0) { // seems slower
   auto norm = norm_induced_2(M);
   if (norm > 1 + 10 * std::numeric_limits<double>::epsilon())
    TRIQS_RUNTIME_ERROR << " Internal Error: norm  >1 !" << norm << r.second;
   n->cache.matrix_lnorms[b] = -std::log(norm);
  }
  }
 return {b3, std::move(M)};
}

//---------------- cache update ------------------------------------

void atomic_correlators_worker::cache_update() {
#ifdef EXT_DEBUG
 tree.graphviz(std::ofstream("tree_start_update_cache"));
#endif
 update_cache_impl(tree.get_root());
}

// -------------------  -------------

void atomic_correlators_worker::update_cache_impl(node n) {

 if ((n == nullptr) || (!n->modified)) return;
 if (n->soft_deleted) TRIQS_RUNTIME_ERROR << " Internal Error: soft deleted node in cache update ";
 update_cache_impl(n->left);
 update_cache_impl(n->right);
 for (int b = 0; b < n_blocks; ++b) { 
  auto r = compute_block_table_and_bound (n, b, -1);
  n->cache.block_table[b] = r.first;
  n->cache.matrix_lnorms[b] = r.second;
  n->cache.matrix_norm_are_valid[b] = false;
 } // n->modified = false;
}

//----------------------------------------------------

atomic_correlators_worker::trace_t atomic_correlators_worker::compute_trace(double epsilon, bool estimator_only) {

 if (tree_size == 0) return sosp->partition_function(config->beta()); // simplifies later code

#ifdef EXT_DEBUG
 tree.graphviz(std::ofstream("tree_start_compute_trace"));
#endif
 auto root = tree.get_root();
 double dt = config->beta() - tree.min_key() + tree.max_key(); // beta - tmax + tmin ! the tree is in REVERSE order

#ifdef EXT_DEBUG
 std::cout << " Trace compu ---------------" << std::endl;
 tree.print(std::cout);
 std::cout << "dt = " << dt << std::endl;
 std::cout << *config << std::endl;
#endif

 auto log_epsilon = -std::log(epsilon);
 std::vector<std::pair<double, int>> to_sort; //(n_blocks);
 double lnorm_threshold = std::numeric_limits<double>::max() - 100;

 for (int b = 0; b < n_blocks; ++b) {
  auto block_lnorm_pair = compute_block_table_and_bound(root, b, lnorm_threshold);
  if (block_lnorm_pair.first == b) { // final structural check B ---> returns to B.
   double l_norm = block_lnorm_pair.second + dt * get_block_emin(b);
   lnorm_threshold = std::min(lnorm_threshold, l_norm + log_epsilon);
   to_sort.emplace_back(l_norm, b);
  }
 }
 // Now sort the blocks non structurally 0 according to the bound
 std::sort(to_sort.begin(), to_sort.end());

 if (estimator_only) {
  double esti = 0;
  for (auto const& x : to_sort) esti += std::exp(-x.first);
  return esti;
 }

 // loop on block, according to estimator, and truncation as Trace_epsilon
 trace_t full_trace = 0;

 for (auto const& b_b : to_sort) { // e_b is a tuple (bound, block_number)

  if (use_truncation && (std::exp(-b_b.first) <= std::abs(full_trace) * epsilon)) break;
  int block_index = b_b.second;

  // computes the matrices, recursively along the modified path in the tree
  auto b_mat = compute_matrix(root, block_index);
  if (b_mat.first == -1) TRIQS_RUNTIME_ERROR << " Internal error : B = -1 after compute matrix : " << block_index;

#ifdef CHECK_AGAINST_LINEAR_COMPUTATION
  auto mat2 = check_one_block_matrix_linear(root, block_index, false);
  if (max_element(abs(b_mat.second - mat2)) > 1.e-10) TRIQS_RUNTIME_ERROR << " Matrix failed against linear computation";
#endif

  // trace (mat * exp(- H * (beta - tmax)) * exp (- H * tmin)) to handle the piece outside of the first-last operators.
  auto dim = get_block_dim(block_index);
  trace_t trace_partial = 0;
  for (int u = 0; u < dim; ++u) trace_partial += b_mat.second(u, u) * std::exp(-dt * get_block_eigenval(block_index, u));

  full_trace += trace_partial; // sum for all blocks
 }

 return full_trace;
}

// code for check/debug
#include "./atomic_correlators_worker.checks.cpp"

} // namespace
