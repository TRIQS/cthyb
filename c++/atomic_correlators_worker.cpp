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

double double_max = std::numeric_limits<double>::max(); // easier to read

// --------------- Computation of the matrix norm --> move into triqs::arrays ------------------------

// this norm is too slow, need to change to another norm...
// double norm_induced_2_impl(matrix_view<double> A, matrix_view<double> B) {
double norm_induced_2_impl(matrix<double> const& A, matrix<double> const& B) {
 // WORKAROUND BUG !!
 // std::cout << A << B << std::endl;
 auto M = A * B;
 triqs::arrays::linalg::eigenelements_worker<matrix_view<double>, true> w(M());
 w.invoke();
 auto const& Es = w.values();
 return std::sqrt(Es(first_dim(Es) - 1)); // ordered is guaranteed by lapack
}

double norm_induced_2(matrix<double> const& A) {
 return (first_dim(A) < second_dim(A) ? norm_induced_2_impl(A, A.transpose()) : norm_induced_2_impl(A.transpose(), A));
}

// -----------------------------------------------

namespace cthyb_matrix {

atomic_correlators_worker::atomic_correlators_worker(configuration& c, sorted_spaces const& sosp_, utility::parameters const& p)
   : config(&c), sosp(&sosp_) {

 // Taking parameters from the inputs
 use_truncation = p["use_truncation"];

 std::string ms = p["trace_estimator"];
 try {
  method =
      std::map<std::string, method_t>{{"FullTrace", method_t::FullTrace}, 
                                      {"EstimateTruncEps", method_t::EstimateTruncEps}}.at(ms);
 }
 catch (...) {
  TRIQS_RUNTIME_ERROR << "Trace method " << ms << " not recognized.";
 }

 make_histograms = p["make_path_histograms"];
 if (make_histograms) {
  for (int i = 0; i < n_orbitals; ++i) histo_opcount.emplace_back(100, "histo_opcount" + std::to_string(i) + ".dat");
 }
}

//------------------------------------------------------------------------------

atomic_correlators_worker::trace_t atomic_correlators_worker::estimate() {
 if (method == method_t::FullTrace) return compute_trace(1.e-15);
 //if (method == method_t::EstimateTruncEps) 
 return compute_trace(0.333); // using epsilon = 1/3 for quick estimator
}

//------------------------------------------------------------------------------

atomic_correlators_worker::trace_t atomic_correlators_worker::full_trace_over_estimator() {
 trace_t r = 1;
 if (method == method_t::EstimateTruncEps) {
  trace_t ft = compute_trace(1.e-15);
  trace_t est = estimate();
  r = ft / est;
  if (!std::isfinite(r)) TRIQS_RUNTIME_ERROR << " full_trace_over_estimator : r not finite" << r << " " << ft << " " << est;
 }
 if (make_histograms) histo_trace_over_estimator << r;
 return r;
}

// ------------------- computation block table (only) -------------

// for subtree at node n, returns B'
// precondition: b !=-1, n != null
// returns -1 if cancellation structural
int atomic_correlators_worker::compute_block_table(node n, int b) {

 if (b < 0) TRIQS_RUNTIME_ERROR << " b <0";
 if (!n->modified) return n->cache.block_table[b];

 int b1 = (n->right ? compute_block_table(n->right, b) : b);
 if (b1 < 0) return b1;

 int b2 = (n->soft_deleted ? b1 : get_op_block_map(n, b1));
 if (b2 < 0) return b2;

 return (n->left ? compute_block_table(n->left, b2) : b2);
}
// ------------------- computation block table and bounds -------------

// for subtree at node n, return (B', bound)
// precondition: b !=-1, n != null
// returns -1 if cancellation structural, -2 if due to threshold
std::pair<int, double> atomic_correlators_worker::compute_block_table_and_bound(node n, int b, double lnorm_threshold) {

 if (b < 0) TRIQS_RUNTIME_ERROR << " b <0";
 // if (n == nullptr) TRIQS_RUNTIME_ERROR << " null ptr";
 if (!n->modified) return {n->cache.block_table[b], n->cache.matrix_lnorms[b]};

 double lnorm = 0;

 int b1 = b;
 if (n->right) {
  std::tie(b1, lnorm) = compute_block_table_and_bound(n->right, b, lnorm_threshold);
  if (b1 < 0) return {b1, 0};
  lnorm += n->cache.dtr * get_block_emin(b1);
 }
 if (lnorm > lnorm_threshold) return {-2, 0};

 int b2 = (n->soft_deleted ? b1 : get_op_block_map(n, b1));
 if (b2 < 0) return {b2, 0};

 int b3 = b2;
 if (n->left) {
  lnorm += n->cache.dtl * get_block_emin(b2);
  if (lnorm > lnorm_threshold) return {-2, 0};
  double lnorm3;
  std::tie(b3, lnorm3) = compute_block_table_and_bound(n->left, b2, lnorm_threshold);
  if (b3 < 0) return {b3, 0};
  lnorm += lnorm3;
 }

 if (lnorm > lnorm_threshold) return {-2, 0};
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

 if (n->right) { // M <- M * exp * r[b]
  dtau_r = double(n->key - tree.min_key(n->right));
  auto dim = second_dim(M); // same as get_block_dim(b2);
  for (int i = 0; i < dim; ++i) M(_, i) *= std::exp(-dtau_r * get_block_eigenval(b1, i));
  if ((first_dim(r.second) == 1) && (second_dim(r.second) == 1))
   M *= r.second(0, 0);
  else
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
  if ((first_dim(l.second) == 1) && (second_dim(l.second) == 1))
   M *= l.second(0, 0);
  else
   M = l.second * M;
 }

 if (updating) {
  n->cache.matrices[b] = M;
  n->cache.matrix_norm_are_valid[b] = true;

  // improve the norm
  if (use_norm_of_matrices_in_cache) { // seems slower
   auto norm = norm_induced_2(M);
   if (norm > 1.000000001) TRIQS_RUNTIME_ERROR << " Internal Error: norm  >1 !" << norm << M;
   n->cache.matrix_lnorms[b] = -std::log(norm);
  }
 }

 return {b3, std::move(M)};
}

//---------------- cache update ------------------------------------

void atomic_correlators_worker::update_cache() {
 update_cache_impl(tree.get_root());
 // analysis
 if (make_histograms) {
  histo_opcount_total << config->size() / 2;
  std::vector<int> opcount(n_orbitals, 0); // maximum number of orbitals is n_orbitals
  for (auto const& p : *config) opcount[p.second.linear_index]++;
  for (int i = 0; i < n_orbitals; ++i) histo_opcount[i] << opcount[i] / 2;
 }
}

// --------------------------------

void atomic_correlators_worker::update_cache_impl(node n) {

 if ((n == nullptr) || (!n->modified)) return;
 if (n->soft_deleted) TRIQS_RUNTIME_ERROR << " Internal Error: soft deleted node in cache update ";
 update_cache_impl(n->left);
 update_cache_impl(n->right);
 n->cache.dtr = (n->right ? double(n->key - tree.min_key(n->right)) : 0);
 n->cache.dtl = (n->left ? double(tree.max_key(n->left) - n->key) : 0);
 for (int b = 0; b < n_blocks; ++b) {
  auto r = compute_block_table_and_bound(n, b, double_max); 
  n->cache.block_table[b] = r.first;
  n->cache.matrix_lnorms[b] = r.second;
  n->cache.matrix_norm_are_valid[b] = false;
 }
 // n->modified = false;
}

// --------------------------------
void atomic_correlators_worker::update_dt(node n) {
 if ((n == nullptr) || (!n->modified)) return;
 update_dt(n->left);
 update_dt(n->right);
 n->cache.dtr = (n->right ? double(n->key - tree.min_key(n->right)) : 0);
 n->cache.dtl = (n->left ? double(tree.max_key(n->left) - n->key) : 0);
}

//----------------------------------------------------

atomic_correlators_worker::trace_t atomic_correlators_worker::compute_trace(double epsilon) {

 if (tree_size == 0) return sosp->partition_function(config->beta()); // simplifies later code

 auto root = tree.get_root();
 double dt = config->beta() - tree.min_key() + tree.max_key(); // beta - tmax + tmin ! the tree is in REVERSE order

#ifdef EXT_DEBUG
 std::cout << " Trace compu ---------------" << std::endl;
 tree.print(std::cout);
 std::cout << "dt = " << dt << std::endl;
 std::cout << *config << std::endl;
 tree.graphviz(std::ofstream("tree_start_compute_trace"));
#endif

 auto log_epsilon = -std::log(epsilon);
 std::vector<std::pair<double, int>> to_sort1, to_sort;
 double lnorm_threshold = double_max - 100;

 update_dt(root); // recompute the dt for modified nodes

 for (int b = 0; b < n_blocks; ++b) {
  auto block_lnorm_pair = compute_block_table_and_bound(root, b, lnorm_threshold);

  if (block_lnorm_pair.first == b) { // final structural check B ---> returns to B.
   double lnorm = block_lnorm_pair.second + dt * get_block_emin(b);
   lnorm_threshold = std::min(lnorm_threshold, lnorm + log_epsilon);
   to_sort1.emplace_back(lnorm, b);
  }
 }

 // recut since lnorm_threshold has evolved during the previous computation
 for (auto const& b_b : to_sort1)
  if (b_b.first <= lnorm_threshold) to_sort.push_back(b_b);

 if (make_histograms) {
  histo_nblock_at_root << to_sort.size();
 }

 if (to_sort.size() == 0) return 0.0; // structural 0

 // Now sort the blocks non structurally 0 according to the bound
 std::sort(to_sort.begin(), to_sort.end());

 // loop on block, according to estimator, and truncation as Trace_epsilon
 trace_t full_trace = 0, first_term = 0;
 auto trace_contrib_block = std::vector<std::pair<double, int>>{};

 int n_bl = to_sort.size();                        // number of blocks
 auto bound_cumul = std::vector<double>(n_bl + 1); // cumulative sum of the bounds
 bound_cumul[n_bl] = 0;
 for (int i = n_bl - 1; i >= 0; --i)
 bound_cumul[i] = bound_cumul[i + 1] + std::exp(-to_sort[i].first) * get_block_dim(to_sort[i].second);

 int bl;
 for (bl = 0; bl < n_bl; ++bl) { // sum over all blocks

  // stopping criterion 
  if (use_truncation && (bl > 0) && (bound_cumul[bl] <= std::abs(full_trace) * epsilon)) break;

  int block_index = to_sort[bl].second;

  // computes the matrices, recursively along the modified path in the tree
  auto b_mat = compute_matrix(root, block_index);
  if (b_mat.first == -1) TRIQS_RUNTIME_ERROR << " Internal error : B = -1 after compute matrix : " << block_index;

#ifdef CHECK_AGAINST_LINEAR_COMPUTATION
  auto mat2 = check_one_block_matrix_linear(root, block_index, false);
  if (max_element(abs(b_mat.second - mat2)) > 1.e-10) TRIQS_RUNTIME_ERROR << " Matrix failed against linear computation";
#endif

  // trace (mat * exp(- H * (beta - tmax)) * exp (- H * tmin)) to handle the piece outside of the first-last operators.
  trace_t trace_partial = 0;
  auto dim = get_block_dim(block_index);
  for (int u = 0; u < dim; ++u) trace_partial += b_mat.second(u, u) * std::exp(-dt * get_block_eigenval(block_index, u));

  if (std::abs(trace_partial) > 1.000001 * dim * std::exp(-to_sort[bl].first))
   TRIQS_RUNTIME_ERROR << "Matrix not bounded by the bound !!!";

  full_trace += trace_partial; // sum for all blocks

  // Analysis
  if (make_histograms) {
   histo_trace_over_bound << std::abs(trace_partial) / std::exp(-to_sort[bl].first);
   trace_contrib_block.emplace_back(std::abs(trace_partial), block_index);
   if (bl == 1) {
    first_term = trace_partial;
    histo_dominant_block_bound << block_index;
    histo_dominant_block_energy_bound << get_block_emin(block_index);
   } else
    histo_trace_first_over_sec_term << trace_partial / first_term;
  }
 } // loop on block

  // Analysis
 if (make_histograms) {
  std::sort(trace_contrib_block.begin(), trace_contrib_block.end(), std::c14::greater<>());
  histo_dominant_block_trace << begin(trace_contrib_block)->second;
  histo_dominant_block_energy_trace << get_block_emin(begin(trace_contrib_block)->second);
  histo_n_block_kept << bl;
  histo_trace_first_term_trace << std::abs(first_term) / std::abs(full_trace);
 }

 if (!std::isfinite(full_trace)) TRIQS_RUNTIME_ERROR << " full_trace not finite" << full_trace;

 return full_trace;
}

// code for check/debug
#include "./atomic_correlators_worker.checks.cpp"

} // namespace
