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

// -----------------------------------------------

namespace cthyb {

// -------- Constructor --------
impurity_trace::impurity_trace(configuration& c, atom_diag const& h_diag_, solve_parameters_t const& p)
   : config(&c),
     h_diag(&h_diag_),
     histo(p.performance_analysis ? new histograms_t(h_diag_.n_blocks()) : nullptr),
     atomic_z(partition_function(*h_diag, config->beta())),
     atomic_norm(0),
     atomic_rho(n_blocks),
     density_matrix(n_blocks) {

 use_norm_as_weight = p.use_norm_as_weight;
 measure_density_matrix = p.measure_density_matrix;
 // init density_matrix block + bool
 for (int bl = 0; bl < n_blocks; ++bl) density_matrix[bl] = bool_and_matrix{false, matrix<double>(get_block_dim(bl), get_block_dim(bl))};

 // prepare atomic_rho and atomic_norm
 if(use_norm_as_weight) {
  auto rho = atomic_density_matrix(h_diag_, config->beta());
  for (int bl = 0; bl < n_blocks; ++bl) {
   atomic_rho[bl] = bool_and_matrix{true, rho[bl] * atomic_z};
   for(int u = 0; u < get_block_dim(bl); ++u) {
    auto xx = rho[bl](u, u);
    atomic_norm += xx * xx;
   }
  }
  atomic_norm = std::sqrt(atomic_norm);
 }
}

//====== Recursive operations ======

// For all recursive operations, the cache on the current node is updated as follows:
//
//               node current: b_current
//               /          \
//              /            \
//             /              \
//            /                \
//           node: b_left     node right: b_right
// tau=beta  <========<=========<   tau=0
// That is, quantities are always updated from right -> left,
// following the time-ordering of beta (left) <- 0 (right).
// Accordingly, the blocks are connected as follows: b_left <- b_current <- b_right

// ------- Computation of the block table (only) -------------

// for subtree at node n, returns B' that block b at the node closest to tau=0 connects to
// precondition: b != -1, n != null
// returns -1 if cancellation structural
int impurity_trace::compute_block_table(node n, int b) {

 if (b < 0) TRIQS_RUNTIME_ERROR << " b < 0";
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

 if (b < 0) TRIQS_RUNTIME_ERROR << " b < 0";
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

// returns {block that b connects to at this node, matrix for this block on node n (if not structurally zero, i.e. if B' != -1)}
std::pair<int, arrays::matrix<double>> impurity_trace::compute_matrix(node n, int b) {

 if (b == -1) return {-1, {}};
 if (n == nullptr) return {b, {}};
 if (!n->modified && n->cache.matrix_norm_valid[b]) return {n->cache.block_table[b], n->cache.matrices[b]};
 bool updating = (!n->modified && !n->cache.matrix_norm_valid[b]);

 double dtau_l = 0, dtau_r = 0;
 auto _ = arrays::range();

 auto r = compute_matrix(n->right, b);
 int b1 = r.first; // exit block of right subtree
 if (b1 == -1) return {-1, {}};

 int b2 = (n->delete_flag ? b1 : get_op_block_map(n, b1)); // relevant block on current node
 if (b2 == -1) return {-1, {}};

 matrix<double> M = (!n->delete_flag ? get_op_block_matrix(n, b1) : make_unit_matrix<double>(get_block_dim(b1)));

 if (n->right) { // M <- M * exp * r[b]
  dtau_r = double(n->key - tree.min_key(n->right));
  auto dim = second_dim(M); // same as get_block_dim(b2);
  for (int i = 0; i < dim; ++i) M(_, i) *= std::exp(-dtau_r * get_block_eigenval(b1, i)); // Create time-evolution matrix e^-H(t'-t)
  if ((first_dim(r.second) == 1) && (second_dim(r.second) == 1))
   M *= r.second(0, 0);
  else
   M = M * r.second; // FIXME could try to optimise lapack call?
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
  n->cache.matrix_norm_valid[b] = true;

  // improve the norm if calculating the full_trace
  if (use_norm_of_matrices_in_cache) { // seems slower
   auto norm = frobenius_norm(M);
   n->cache.matrix_lnorms[b] = -std::log(norm);
   if (!std::isfinite(-std::log(norm))) {
    n->cache.matrix_lnorms[b] = double_max;
   }
  }
 }

 return {b3, std::move(M)};
}

// ------- Update the cache -----------------------

void impurity_trace::update_cache() {
 update_cache_impl(tree.get_root());
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
 // This is not necessary here as all modified nodes are "cleared"
 //  by tree::clear_modified in the try/cancel/confirm set
 // n->modified = false;
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
// Returns MC atomic weight and reweighting = trace/(atomic weight)
std::pair<double, impurity_trace::trace_t> impurity_trace::compute(double p_yee, double u_yee) {

 double epsilon = 1.e-15; // Machine precision
 auto log_epsilon0 = -std::log(1.e-15);
 double lnorm_threshold = double_max - 100;
 std::vector<std::pair<double, int>> init_to_sort_lnorm_b, to_sort_lnorm_b; // pairs of lnorm and b to sort in order of bound

 // simplifies later code
 if (tree_size == 0) {
  if(use_norm_as_weight) {
   density_matrix = atomic_rho;
   return {atomic_norm, atomic_z / atomic_norm};
  } else return {atomic_z, 1};
 }

 auto root = tree.get_root();
 // beta - tmax + tmin ! the tree is in REVERSE order
 double dtau_beta = config->beta() - tree.min_key();
 double dtau_0 = double(tree.max_key());
 double dtau = dtau_beta + dtau_0;

//FIXME
// #ifdef EXT_DEBUG
//  std::cout << " Trace computed ---------------" << std::endl;
//  tree.print(std::cout);
//  std::cout << "dtau = " << dtau << std::endl;
//  std::cout << *config << std::endl;
//  tree.graphviz(std::ofstream("tree_start_compute_trace"));
// #endif

 update_dtau(root); // recompute the dtau for modified nodes

 for (int b = 0; b < n_blocks; ++b) {
  auto block_lnorm_pair = compute_block_table_and_bound(root, b, lnorm_threshold);

  // Check that the final block is the same as the initial block or -1, indicating structural cancellation
  // This guarantees that the density matrix is blockwise diagonal (otherwise the code will have thrown an error).
  if (measure_density_matrix) {
   if ((block_lnorm_pair.first != b) && (block_lnorm_pair.first != -1))
    TRIQS_RUNTIME_ERROR << "The product of atomic operators has a matrix element in the off-diagonal block ("
                        << b << "," << block_lnorm_pair.first << ")\n" << *config;
  }

  if (block_lnorm_pair.first == b) { // final structural check B ---> returns to B.
   double lnorm = block_lnorm_pair.second + dtau * get_block_emin(b);
   lnorm_threshold = std::min(lnorm_threshold, lnorm + log_epsilon0);
   init_to_sort_lnorm_b.emplace_back(lnorm, b);
  }
 }

 // recut since lnorm_threshold evolved in the previous loop
 for (auto const& b_b : init_to_sort_lnorm_b)
  if (b_b.first <= lnorm_threshold) to_sort_lnorm_b.push_back(b_b);

 if (histo) histo->n_block_at_root << to_sort_lnorm_b.size();

 if (to_sort_lnorm_b.size() == 0) return {0.0, 1}; // structural 0

 // Now sort the blocks non structurally 0 according to the bound
 std::sort(to_sort_lnorm_b.begin(), to_sort_lnorm_b.end());

 // Prepare to loop over all blocks (in sorted order).
 // According to estimator, truncate as epsilon.
 trace_t full_trace = 0, first_term = 0;
 double norm_trace_sq = 0, trace_abs =0;

 // Put density_matrix to "not recomputed"
 for (int bl = 0; bl < n_blocks; ++bl) density_matrix[bl].is_valid = false;

 auto trace_contrib_block = std::vector<std::pair<double, int>>{};

 int n_bl = to_sort_lnorm_b.size();                // number of blocks
 auto bound_cumul = std::vector<double>(n_bl + 1); // cumulative sum of the bounds
 // The contribution to the trace from block B is bounded: |Tr_B| <= dim(B) * sum_{B} e^{Emin(B)*dtau}
 // Here we calculate the cumulative bound from each contributing (structurally non-zero) block to
 // determine at which block we have exceeded the bound and hence can stop.
 // Can tighten bound on trace by using sqrt(dim(B)) in the case of Frobenius norm only.
 bound_cumul[n_bl] = 0;
 if (!use_norm_as_weight) {
  for (int bl = n_bl - 1; bl >= 0; --bl)
   bound_cumul[bl] = bound_cumul[bl + 1] + std::exp(-to_sort_lnorm_b[bl].first) * std::sqrt(get_block_dim(to_sort_lnorm_b[bl].second));
 } else {
  for (int bl = n_bl - 1; bl >= 0; --bl) bound_cumul[bl] = bound_cumul[bl + 1] + std::exp(-to_sort_lnorm_b[bl].first);
 }

 int bl;
 for (bl = 0; bl < n_bl; ++bl) { // sum over all blocks

  // stopping criterion
  if ((bl > 0) && (bound_cumul[bl] <= std::abs(full_trace) * epsilon)) break;

  int block_index = to_sort_lnorm_b[bl].second; // index in original (unsorted) order

  // additionnal Yee quick return criterion
  if (p_yee >= 0.0) {
   auto current_weight = (use_norm_as_weight ? std::sqrt(norm_trace_sq) : full_trace);
   auto pmax = std::abs(p_yee) * (std::abs(current_weight) + bound_cumul[bl]);
   if (pmax < u_yee) return {0, 1}; // pmax < u, we can reject
  }

  // computes the matrices, recursively along the modified path in the tree
  auto b_mat = compute_matrix(root, block_index); // b_mat = {block that b connects to, matrix for this block}
  if (b_mat.first == -1) TRIQS_RUNTIME_ERROR << " Internal error : B = -1 after compute matrix : " << block_index;

#ifdef CHECK_AGAINST_LINEAR_COMPUTATION
  auto b_mat2 = check_one_block_matrix_linear(root, block_index, false);
  if (max_element(abs(b_mat.second - b_mat2)) > 1.e-10) TRIQS_RUNTIME_ERROR << " Matrix failed against linear computation";
#endif

  // trace(mat * exp(- H * (beta - tmax)) * exp (- H * tmin)) to handle the piece outside of the first-last operators.
  trace_t trace_partial = 0;
  auto dim = get_block_dim(block_index);
  for (int u = 0; u < dim; ++u) {
   auto x = b_mat.second(u, u) * std::exp(-dtau * get_block_eigenval(block_index, u));
   trace_partial += x;
   trace_abs += std::abs(x);
  }

  if (use_norm_as_weight) { // else we are not allowed to compute this matrix, may make no sense
   // recompute the density matrix
   density_matrix[block_index].is_valid = true;
   double norm_trace_sq_partial = 0;
   auto& mat = density_matrix[block_index].mat;
   for (int u = 0; u < dim; ++u) {
    for (int v = 0; v < dim; ++v) {
     mat(u, v) = b_mat.second(u, v) *
                 std::exp(-dtau_beta * get_block_eigenval(block_index, u) - dtau_0 * get_block_eigenval(block_index, v));
     double xx = std::abs(mat(u, v));
     norm_trace_sq_partial += xx * xx;
    }
   }
   norm_trace_sq += norm_trace_sq_partial;
   // internal check
   if (std::abs(trace_partial) - 1.0000001 * std::sqrt(norm_trace_sq_partial) * get_block_dim(block_index) > 1.e-15)
    TRIQS_RUNTIME_ERROR << "|trace| > dim * norm" << trace_partial << " " << std::sqrt(norm_trace_sq_partial) << "  " << trace_abs;
   if (std::abs(trace_partial - trace(mat)) > 1.e-15) TRIQS_RUNTIME_ERROR << "Internal error : trace and density mismatch";
  }

#ifdef CHECK_MATRIX_BOUNDED_BY_BOUND
  if (std::abs(trace_partial) > 1.000001 * dim * std::exp(-to_sort_lnorm_b[bl].first))
   TRIQS_RUNTIME_ERROR << "Matrix not bounded by the bound ! test is " << std::abs(trace_partial) <<" < " << dim * std::exp(-to_sort_lnorm_b[bl].first);
#endif

  full_trace += trace_partial; // sum for all blocks

  // Analysis
  if (histo) {
   histo->trace_over_bound << std::abs(trace_partial) / std::exp(-to_sort_lnorm_b[bl].first);
   trace_contrib_block.emplace_back(std::abs(trace_partial), block_index);
   if (bl == 1) {
    first_term = trace_partial;
    histo->dominant_block_bound << block_index;
    histo->dominant_block_energy_bound << get_block_emin(block_index);
   } else
    histo->trace_first_over_sec_term << trace_partial / first_term;
  }
 } // loop on block

 double norm_trace = std::sqrt(norm_trace_sq);
 if (!std::isfinite(full_trace)) TRIQS_RUNTIME_ERROR << " full_trace not finite" << full_trace;

 // Analysis
 if (histo) {
  histo->trace_over_norm << std::abs(full_trace)/ norm_trace;
  histo->trace_abs_over_norm << trace_abs/ norm_trace;
  histo->trace_over_trace_abs << full_trace / trace_abs;
  std::sort(trace_contrib_block.begin(), trace_contrib_block.end(), std::c14::greater<>());
  histo->dominant_block_trace << begin(trace_contrib_block)->second;
  histo->dominant_block_energy_trace << get_block_emin(begin(trace_contrib_block)->second);
  histo->n_block_kept << bl;
  histo->trace_first_term_trace << std::abs(first_term) / std::abs(full_trace);
 }

 // return {weight, reweighting}
 if (!use_norm_as_weight) return {full_trace, 1};
 // else determine reweighting
 auto rw = full_trace / norm_trace;
 if (!std::isfinite(rw)) rw = 1;
 //FIXME if (!std::isfinite(rw)) TRIQS_RUNTIME_ERROR << "Atomic correlators : reweight not finite" << full_trace << " "<< norm_trace;
 return {norm_trace, rw};
 }

// code for check/debug
#include "./impurity_trace.checks.cpp"

} // namespace
