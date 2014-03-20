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

//double norm_induced_2_impl(matrix_view<double> A, matrix_view<double> B) {
double norm_induced_2_impl(matrix<double> const &A, matrix<double> const & B) {
 // WORKAROUND BUG !!
 //std::cout << A << B << std::endl;
 auto M = A * B;
 triqs::arrays::linalg::eigenelements_worker<matrix_view<double>, true> w(M());
 w.invoke();
 auto const& Es = w.values();
 return std::sqrt(Es(first_dim(Es) - 1)); // ordered is guaranteed by lapack
}

double norm_induced_2(matrix<double> const& A) {
 return (first_dim(A) < second_dim(A) ? norm_induced_2_impl(A, A.transpose()) : norm_induced_2_impl(A.transpose(), A));
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

  //for (int i = 0; i < 20; ++i) histo_n_blocks_after_steps.emplace_back(sosp->n_subspaces(), "histo_n_blocks_after_steps" + std::to_string(i) + ".dat");
 
  for (int i = 0; i < n_orbitals; ++i) 
   histo_opcount.emplace_back(100, "histo_opcount" + std::to_string(i) +  ".dat");

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
 //if (make_histograms) histos["FullTrace_over_Estimate"] << std::abs(r);
 return r;
}

// ------------------- computation block table and bounds -------------

// for subtree at node n, return (B', bound)
// precondition: b !=-1, n != null
// returns -1 if cancellation structural, -2 if due to threshold
int atomic_correlators_worker::compute_block_table(node n, int b) {

 if (b < 0) TRIQS_RUNTIME_ERROR << " b <0";
 if (!n->modified) return n->cache.block_table[b];

 int b1 = b;
 if (n->right) {
  b1 = compute_block_table(n->right, b);
  if (b1 < 0) return b1;
 }

 int b2 = (n->soft_deleted ? b1 : get_op_block_map(n, b1));
 if (b2 < 0) return b2;

 int b3 = b2;
 if (n->left) b3 = compute_block_table(n->left, b2);
 return b3;
}
// ------------------- computation block table and bounds -------------

// for subtree at node n, return (B', bound)
// precondition: b !=-1, n != null
// returns -1 if cancellation structural, -2 if due to threshold
std::pair<int, double> atomic_correlators_worker::compute_block_table_and_bound(node n, int b, double lnorm_threshold) { 

 if (b < 0) TRIQS_RUNTIME_ERROR << " b <0"; 
 //if (n == nullptr) TRIQS_RUNTIME_ERROR << " null ptr";
 if (!n->modified) return {n->cache.block_table[b], n->cache.matrix_lnorms[b]};
 
 histo_appel << b;
 n_call++;

 double lnorm= 0;

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

 //auto ln1 = norm_induced_2(M);

 if (n->right) { // M <- M* exp * r[b]
  dtau_r = double(n->key - tree.min_key(n->right));
  auto dim = second_dim(M); // same as get_block_dim(b2);
  for (int i = 0; i < dim; ++i) M(_, i) *= std::exp(-dtau_r * get_block_eigenval(b1, i));
  //histo_dims_1x1 << ((first_dim(r.second) == 1) && (second_dim(r.second) == 1));
  if (make_histograms) histo_dims_1x1 << ((first_dim(M) == 1) && (second_dim(M) == 1));
  if ((first_dim(r.second) == 1) && (second_dim(r.second) == 1))
   M *= r.second(0, 0);
  else
   M = M * r.second; // optimise lapack call ?
  //auto ratio =  norm_induced_2(M) /( ln1 * norm_induced_2(r.second));
  //if (ratio < 0.01) {
   //std::cout << ratio << " " << norm_induced_2(M) << "  " << ln1 * norm_induced_2(r.second) << std::endl;
   //std::cout << M << r.second << std::endl;
  //}
  //histo_perte_prod << ratio;
 }


 int b3 = b2;
 if (n->left) { // M <- l[b] * exp * M
  auto l = compute_matrix(n->left, b2);
  b3 = l.first;
  if (b3 == -1) return {-1, {}};
  dtau_l = double(tree.max_key(n->left) - n->key);
  auto dim = first_dim(M); // same as get_block_dim(b1);
  for (int i = 0; i < dim; ++i) M(i, _) *= std::exp(-dtau_l * get_block_eigenval(b2, i));
  //histo_dims_1x1 << ((first_dim(l.second) == 1) && (second_dim(l.second) == 1));
  if (make_histograms) histo_dims_1x1 << ((first_dim(M) == 1) && (second_dim(M) == 1));
  if ((first_dim(l.second) == 1) && (second_dim(l.second) == 1))
   M *= l.second(0, 0);
  else
   M = l.second * M;
  //histo_perte_prod << norm_induced_2(M) / (ln1 * norm_induced_2(l.second));
 }

 if (updating) {
  n->cache.matrices[b] = M;
  n->cache.matrix_norm_are_valid[b] = true;

  // improve the norm
  if (use_norm_of_matrices_in_cache) { // seems slower
   // std::cout  << first_dim(M) << " "<< second_dim(M) << std::endl;
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

 // analysis
 if (make_histograms) {
  histo_opcount_total << config->size()/2;
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
  auto r = compute_block_table_and_bound(n, b, double_max);//, 0);
  n->cache.block_table[b] = r.first;
  n->cache.matrix_lnorms[b] = r.second;
  n->cache.matrix_norm_are_valid[b] = false;
  //if (b>10) break;
 } // n->modified = false;
}

// --------------------------------
//// ---- TRY update dt
void atomic_correlators_worker::update_dt(node n) {
 if ((n == nullptr) || (!n->modified)) return;
 update_dt(n->left);
 update_dt(n->right);
 n->cache.dtr = (n->right ? double(n->key - tree.min_key(n->right)) : 0);
 n->cache.dtl = (n->left ? double(tree.max_key(n->left) - n->key) : 0);
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
 std::vector<std::pair<double, int>> to_sort1, to_sort; 
 double lnorm_threshold = double_max - 100;

 update_dt(root); // recompute the dt for modified nodes

//#define NO_BOUND
#ifndef NO_BOUND
 for (int b = 0; b < n_blocks; ++b) {
  n_call = 0;
  block_dominant = 0;
  int bl = (block_dominant + b) % n_blocks;
  auto block_lnorm_pair = compute_block_table_and_bound(root, bl, lnorm_threshold);

  if (make_histograms) {
   if (block_lnorm_pair.first == -1) histo_cut1 << b;
   if (block_lnorm_pair.first == -2) histo_cut2 << b;
   if (block_lnorm_pair.first == -2) histo_cut2_energy << get_block_emin(b);
   if (b < 5) histo_appel2 << n_call;
   if (abs(300 - b) < 5) histo_appel3 << n_call;
  }

  if (block_lnorm_pair.first == bl) { // final structural check B ---> returns to B.
   double lnorm = block_lnorm_pair.second + dt * get_block_emin(b);
   //lnorm_threshold = std::min(lnorm_threshold, lnorm + log_epsilon);
   lnorm_threshold = std::min(lnorm_threshold, lnorm + (use_only_first_term_in_trace ? 0 : log_epsilon));
   to_sort1.emplace_back(lnorm, bl);
  }
 }
#else
 for (int b = 0; b < n_blocks; ++b) {
  n_call = 0;
  block_dominant = 0;
  auto bf = compute_block_table(root, b);
  if (bf == b) { // final structural check B ---> returns to B.
   to_sort1.emplace_back(1, b);
  }
 }
#endif


 // recut since lnorm_threshold has evolved during the previous computation
 for (auto const& b_b : to_sort1)
  if (b_b.first <= lnorm_threshold) to_sort.push_back(b_b);

 if (make_histograms)  { 
  histo_nblock_at_root << to_sort.size();
  histo_n_block_examined_first_pass << to_sort1.size();
 }

 if (to_sort.size() == 0) return 0.0; // structural 0

 // Now sort the blocks non structurally 0 according to the bound
 std::sort(to_sort.begin(), to_sort.end());

 if (estimator_only) {
  double esti = 0;
  for (auto const& x : to_sort) esti += std::exp(-x.first);
  return esti;
 }

 // loop on block, according to estimator, and truncation as Trace_epsilon
 trace_t full_trace = 0, first_term = 0;
 int n_block_kept=0; // count the number of block that the trace will keep at the end (after truncation)
 auto trace_contrib_block = std::vector<std::pair<double, int>>{};

 for (auto const& b_b : to_sort) { // e_b is a tuple (bound, block_number)

  if (use_truncation && (n_block_kept>0) && (std::exp(-b_b.first) <= std::abs(full_trace) * epsilon)) break;
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

  if (std::abs(trace_partial) > 1.000001 * dim * std::exp(-b_b.first)) TRIQS_RUNTIME_ERROR << "Matrix not bounded by the bound !!!";

  if (make_histograms) histo_trace_over_bound << std::abs(trace_partial) / std::exp(-b_b.first);

  full_trace += trace_partial; // sum for all blocks
  // Analysis
  n_block_kept++;
  if (n_block_kept == 1) block_dominant = block_index;

  if (make_histograms) {
   trace_contrib_block.emplace_back(std::abs(trace_partial), block_index);
   if (n_block_kept == 1) {
    first_term = trace_partial;
    //std::cout << "dominant block : " << block_index << std::endl;
    histo_dominant_block_bound << block_index;
    histo_dominant_block_energy_bound << get_block_emin(block_index);
   } else
    histo_trace_first_over_sec_term << trace_partial / first_term;
  }
  if (use_only_first_term_in_trace) break; // only the first term in the trace
 } // loop on block

 if (n_block_kept == 0) TRIQS_RUNTIME_ERROR << "no block kept in trace computation !";
 if (make_histograms) {
  std::sort(trace_contrib_block.begin(), trace_contrib_block.end(), std::c14::greater<>());
  histo_dominant_block_trace << begin(trace_contrib_block)->second;
  histo_dominant_block_energy_trace << get_block_emin(begin(trace_contrib_block)->second);

  histo_n_block_kept << n_block_kept;
  histo_trace_first_term_trace << std::abs(first_term) / std::abs(full_trace);
 }
 return full_trace;
}

// code for check/debug
#include "./atomic_correlators_worker.checks.cpp"

} // namespace
