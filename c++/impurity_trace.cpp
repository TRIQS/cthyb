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
#include "impurity_trace.hpp"
#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>
#include <algorithm>
#include <limits>
#include <triqs/arrays/linalg/eigenelements.hpp>

//#define CHECK_ALL
#ifdef CHECK_ALL
#define CHECK_CACHE
#define CHECK_AGAINST_LINEAR_COMPUTATION
#define CHECK_MATRIX_BOUNDED_BY_BOUND
#endif

double double_max = std::numeric_limits<double>::max(); // easier to read

// --------------- FIXME Computation of the matrix norm --> move into triqs::arrays ------------------------

namespace triqs {
namespace arrays {

 double frobenius_norm(matrix<double> const& a) {
  return std::sqrt(fold([](double r, double x)->double {
   auto ab = std::abs(x);
   return r + ab * ab;
  })(a));
 }
}
}

// -----------------------------------------------

namespace cthyb {

// -------- Constructor --------
impurity_trace::impurity_trace(configuration& c, sorted_spaces const& sosp_, solve_parameters_t const& p)
   : config(&c), sosp(&sosp_), histo(p.make_histograms ? new histograms_t(sosp_.n_subspaces()) : nullptr) {

 // Taking parameters from the inputs
 use_trace_estimator = p.use_trace_estimator;
 if (use_trace_estimator) 
  method = method_t::estimate;
 else
  method = method_t::full_trace;

 if (histo) {
  for (int i = 0; i < n_orbitals; ++i) histo->opcount.emplace_back(100, "histo_opcount" + std::to_string(i) + ".dat");
 }
}

// -------- Calculate the estimate of the trace -------------------------------------------

impurity_trace::trace_t impurity_trace::estimate(double p_yee, double u_yee) {
 if (method == method_t::full_trace) return compute_trace(true, p_yee, u_yee);
 //if (method == method_t::estimate)
 return compute_trace(false, p_yee, 0);
}

// -------- Calculate the ratio of the full trace to the estimator ----------------------------

impurity_trace::trace_t impurity_trace::full_trace_over_estimator() {
 trace_t r = 1;
 if (method == method_t::estimate) {
  trace_t ft = compute_trace(true, -1, 0);
  trace_t est = estimate(-1, 0);
  r = ft / est;
  if (std::abs(r - 2.0 / 3.0) > 2.0 / 3.0 + 0.001) TRIQS_RUNTIME_ERROR << " estimator out of bounds !! " << r;
  if (!std::isfinite(r)) {
    // It might be that both full trace and est are zero. This is still OK and r should be 1.
    if (std::abs(est) < std::numeric_limits<trace_t>::epsilon() && std::abs(ft) < std::numeric_limits<trace_t>::epsilon())
      r = 1;
    else
      TRIQS_RUNTIME_ERROR << "full_trace_over_estimator: r not finite" << r << " " << ft << " " << est << " " << *config;
  }
 }
 if (histo) histo->trace_over_estimator << r;
 return r;
}

// ------- Computation of the block table (only) -------------

// for subtree at node n, returns B'
// precondition: b !=-1, n != null
// returns -1 if cancellation structural
int impurity_trace::compute_block_table(node n, int b) {

 if (b < 0) TRIQS_RUNTIME_ERROR << " b <0";
 if (!n->modified) return n->cache.block_table[b];

 int b1 = (n->right ? compute_block_table(n->right, b) : b);
 if (b1 < 0) return b1;

 int b2 = (n->delete_flag ? b1 : get_op_block_map(n, b1));
 if (b2 < 0) return b2;

 return (n->left ? compute_block_table(n->left, b2) : b2);
}
// -------- Computation of the block table and bounds -------------

// for subtree at node n, return (B', bound)
// precondition: b !=-1, n != null
// returns -1 for structural and/or threshold cancellation
std::pair<int, double> impurity_trace::compute_block_table_and_bound(node n, int b, double lnorm_threshold, bool use_threshold) {

 if (b < 0) TRIQS_RUNTIME_ERROR << " b <0";
 // if (n == nullptr) TRIQS_RUNTIME_ERROR << " null ptr";
 if (!n->modified) return {n->cache.block_table[b], n->cache.matrix_lnorms[b]};

 double lnorm = 0;

 int b1 = b;
 if (n->right) {
  std::tie(b1, lnorm) = compute_block_table_and_bound(n->right, b, lnorm_threshold, use_threshold);
  if (b1 < 0) return {b1, 0};
  lnorm += n->cache.dtau_r * get_block_emin(b1);
 }
 if (use_threshold && (lnorm > lnorm_threshold)) return {-1, 0};

 int b2 = (n->delete_flag ? b1 : get_op_block_map(n, b1));
 if (b2 < 0) return {b2, 0};

 int b3 = b2;
 if (n->left) {
  lnorm += n->cache.dtau_l * get_block_emin(b2);
  if (use_threshold && (lnorm > lnorm_threshold)) return {-1, 0};
  double lnorm3;
  std::tie(b3, lnorm3) = compute_block_table_and_bound(n->left, b2, lnorm_threshold, use_threshold);
  if (b3 < 0) return {b3, 0};
  lnorm += lnorm3;
 }

 if (use_threshold && (lnorm > lnorm_threshold)) return {-1, 0};

 if (std::isinf(lnorm)){
  lnorm = double_max;
  if (lnorm<0) TRIQS_RUNTIME_ERROR << "Negative lnorm in compute_block_table_and_bound!";
 }

 return {b3, lnorm};
}

// -------- Computation of the matrix ------------------------------

std::pair<int, arrays::matrix<double>> impurity_trace::compute_matrix(node n, int b) {

 if (b == -1) return {-1, {}};
 if (n == nullptr) return {b, {}};
 if (!n->modified && n->cache.matrix_norm_valid[b]) return {n->cache.block_table[b], n->cache.matrices[b]};
 bool updating = (!n->modified && !n->cache.matrix_norm_valid[b]);

 double dtau_l = 0, dtau_r = 0;
 auto _ = arrays::range();

 auto r = compute_matrix(n->right, b);
 int b1 = r.first;
 if (b1 == -1) return {-1, {}};

 int b2 = (n->delete_flag ? b1 : get_op_block_map(n, b1));
 if (b2 == -1) return {-1, {}};

 matrix<double> M = (!n->delete_flag ? get_op_block_matrix(n, b1) : make_unit_matrix<double>(get_block_dim(b1)));

 if (n->right) { // M <- M * exp * r[b]
  dtau_r = double(n->key - tree.min_key(n->right));
  auto dim = second_dim(M); // same as get_block_dim(b2);
  for (int i = 0; i < dim; ++i) M(_, i) *= std::exp(-dtau_r * get_block_eigenval(b1, i));
  //for (int i = 1; i < dim; ++i) M(_, i) *= std::exp(-dtau_r * (get_block_eigenval(b1, i) - get_block_eigenval(b1, 0)));
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
  //for (int i = 1; i < dim; ++i) M(i, _) *= std::exp(-dtau_l* ( get_block_eigenval(b2, i) - get_block_eigenval(b2, 0))); FIXME
  if ((first_dim(l.second) == 1) && (second_dim(l.second) == 1))
   M *= l.second(0, 0);
  else
   M = l.second * M;
 }
 
 if (updating) {
  n->cache.matrices[b] = M;
  n->cache.matrix_norm_valid[b] = true;

  // improve the norm
  if ((method == method_t::full_trace) && use_norm_of_matrices_in_cache) { // seems slower
   auto norm = frobenius_norm(M);
   n->cache.matrix_lnorms[b] = -std::log(norm);
   if (!std::isfinite(-std::log(norm))) {
    //std::cerr << "norm is not finite "<< norm <std::endl; FIXME
    n->cache.matrix_lnorms[b] = double_max;
   }
  }
 }

 return {b3, std::move(M)};
}

// ------- Update the cache -----------------------

void impurity_trace::update_cache() {
 update_cache_impl(tree.get_root());
 // analysis
 if (histo) {
  histo->opcount_total << config->size() / 2;
  std::vector<int> opcount(n_orbitals, 0); // maximum number of orbitals is n_orbitals
  for (auto const& p : *config) opcount[p.second.linear_index]++;
  for (int i = 0; i < n_orbitals; ++i) histo->opcount[i] << opcount[i] / 2;
 }
}

// --------------------------------

void impurity_trace::update_cache_impl(node n) {

 if ((n == nullptr) || (!n->modified)) return;
 if (n->delete_flag) TRIQS_RUNTIME_ERROR << " Internal Error: node flagged for deletion in cache update ";
 update_cache_impl(n->left);
 update_cache_impl(n->right);
 n->cache.dtau_r = (n->right ? double(n->key - tree.min_key(n->right)) : 0);
 n->cache.dtau_l = (n->left ? double(tree.max_key(n->left) - n->key) : 0);
 for (int b = 0; b < n_blocks; ++b) {
  auto r = compute_block_table_and_bound(n, b, double_max, false);
  n->cache.block_table[b] = r.first;
  n->cache.matrix_lnorms[b] = r.second;
  n->cache.matrix_norm_valid[b] = false;
 }
 // n->modified = false; FIXME
}

// -------- Calculate the dtau for a given node to its left and right neighbours ----------------
void impurity_trace::update_dtau(node n) {
 if ((n == nullptr) || (!n->modified)) return;
 update_dtau(n->left);
 update_dtau(n->right);
 n->cache.dtau_r = (n->right ? double(n->key - tree.min_key(n->right)) : 0);
 n->cache.dtau_l = (n->left ? double(tree.max_key(n->left) - n->key) : 0);
}

//-------- Compute the full trace ------------------------------------------

impurity_trace::trace_t impurity_trace::compute_trace(bool to_machine_precision, double p_yee, double u_yee) {

 double epsilon = (to_machine_precision ? 1.e-15 : 0.333);

 if (tree_size == 0) return sosp->partition_function(config->beta()); // simplifies later code

 auto root = tree.get_root();
 double dtau = config->beta() - tree.min_key() + tree.max_key(); // beta - tmax + tmin ! the tree is in REVERSE order

#ifdef EXT_DEBUG
 std::cout << " Trace compu ---------------" << std::endl;
 tree.print(std::cout);
 std::cout << "dtau = " << dtau << std::endl;
 std::cout << *config << std::endl;
 tree.graphviz(std::ofstream("tree_start_compute_trace"));
#endif

 auto log_epsilon0 = -std::log(1.e-15);
 std::vector<std::pair<double, int>> to_sort1, to_sort;
 double lnorm_threshold = double_max - 100;

 update_dtau(root); // recompute the dtau for modified nodes

 for (int b = 0; b < n_blocks; ++b) {
  auto block_lnorm_pair = compute_block_table_and_bound(root, b, lnorm_threshold);

  if (block_lnorm_pair.first == b) { // final structural check B ---> returns to B.
   double lnorm = block_lnorm_pair.second + dtau * get_block_emin(b);
   lnorm_threshold = std::min(lnorm_threshold, lnorm + log_epsilon0);
   to_sort1.emplace_back(lnorm, b);
  }
 }

 // recut since lnorm_threshold has evolved during the previous computation
 for (auto const& b_b : to_sort1)
  if (b_b.first <= lnorm_threshold) to_sort.push_back(b_b);

 if (histo) {
  histo->n_block_at_root << to_sort.size();
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
  if ((bl > 0) && (bound_cumul[bl] <= std::abs(full_trace) * epsilon)) break;

  int block_index = to_sort[bl].second;

  // additionnal Yee quick return criterion
  if (p_yee >= 0.0) {
   if (to_machine_precision) {
    auto pmax = std::abs(p_yee) * (std::abs(full_trace) + bound_cumul[bl] * get_block_dim(block_index));
    if (pmax < u_yee) return 0; // pmax < u, we can reject
   } else {
    // correct but then the bound of the estimator may not be fullfilled.
    //auto pmin = std::abs(p_yee) * (std::abs(full_trace) - bound_cumul[bl] * get_block_dim(block_index));
    //if (pmin > 1) return full_trace; // pmin > 1 > u, we will accept and GET THIS ESTIMATION OF THE TRACE WHATEVER u
   }
  }

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
  for (int u = 0; u < dim; ++u) trace_partial += b_mat.second(u, u) * std::exp(-dtau * get_block_eigenval(block_index, u));
  //trace_partial *= std::exp(dtau * get_block_eigenval(block_index, 0) - to_sort[bl].first);

#ifdef CHECK_MATRIX_BOUNDED_BY_BOUND
  if (std::abs(trace_partial) > 1.000001 * dim * std::exp(-to_sort[bl].first))
   TRIQS_RUNTIME_ERROR << "Matrix not bounded by the bound ! test is " << std::abs(trace_partial) <<" < " << dim * std::exp(-to_sort[bl].first);
#endif

  full_trace += trace_partial; // sum for all blocks

  // Analysis
  if (histo) {
   histo->trace_over_bound << std::abs(trace_partial) / std::exp(-to_sort[bl].first);
   trace_contrib_block.emplace_back(std::abs(trace_partial), block_index);
   if (bl == 1) {
    first_term = trace_partial;
    histo->dominant_block_bound << block_index;
    histo->dominant_block_energy_bound << get_block_emin(block_index);
   } else
    histo->trace_first_over_sec_term << trace_partial / first_term;
  }
 } // loop on block

  // Analysis
 if (histo) {
  std::sort(trace_contrib_block.begin(), trace_contrib_block.end(), std::c14::greater<>());
  histo->dominant_block_trace << begin(trace_contrib_block)->second;
  histo->dominant_block_energy_trace << get_block_emin(begin(trace_contrib_block)->second);
  histo->n_block_kept << bl;
  histo->trace_first_term_trace << std::abs(first_term) / std::abs(full_trace);
 }

 if (!std::isfinite(full_trace)) TRIQS_RUNTIME_ERROR << " full_trace not finite" << full_trace;

 return full_trace;
}

// code for check/debug
#include "./impurity_trace.checks.cpp"

} // namespace
