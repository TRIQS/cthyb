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
using namespace triqs_cthyb;

//-------------------- Cache integrity check --------------------------------

void impurity_trace::check_cache_integrity(bool print) {
#ifdef CHECK_CACHE
  static int check_counter = 0;
  ++check_counter;
  if (check_counter % 10 == 0) {
    if (print) std::cout << " ---- Cache integrity check ---- " << std::endl;
    if (print) std::cout << " check_counter = config number = " << check_counter << std::endl;
    if (print) tree.graphviz(std::ofstream("tree_cache_check"));
    foreach_subtree_first(tree, [&](node y) { this->check_cache_integrity_one_node(y, print); });
    if (print) std::cout << " ---- Cache integrity completed ---- " << std::endl;
  }
#endif
}

//--------------------- Compute block table for one subtree, using an ordered traversal of the subtree -------------------

int impurity_trace::check_one_block_table_linear(node n, int b, bool print) {

  int B = b;
  foreach_reverse(tree, n, [&](node y) {
    if (B == -1) return;
    auto BB = B;
    B       = (y->delete_flag ? B : this->get_op_block_map(y, B));
    if (print)
      std::cout << "linear computation : " << y->key << " " << y->op.dagger << " " << y->op.linear_index << " | " << BB << " -> " << B << std::endl;
  });
  return B;
}

//--------------------- Compute block table for one subtree, using an ordered traversal of the subtree -------------------

matrix_t impurity_trace::check_one_block_matrix_linear(node top, int b) {

  node p     = tree.max(top);
  matrix_t M = nda::eye<h_scalar_t>(get_block_dim(b));
  auto _     = arrays::range::all;

  foreach_reverse(tree, top, [&](node n) {
    // multiply by the exponential unless it is the first call, i.e. first operator n==p
    if (n != p) {
      auto dtau = double(n->key - p->key);
      //  M <- exp * M
      auto dim = M.shape()[0];
      for (int i = 0; i < dim; ++i) M(i, _) *= std::exp(-dtau * get_block_eigenval(b, i));
      // M <- Op * M
    }
    // multiply by operator matrix unless it is delete_flag
    if (!n->delete_flag) {
      int bp = this->get_op_block_map(n, b);
      if (bp == -1) TRIQS_RUNTIME_ERROR << " Nasty error ";
      M = get_op_block_matrix(n, b) * M;
      b = bp;
    }
    p = n;
  });

  return M;
}
//-------------------- Cache integrity check for one node --------------------------------

void impurity_trace::check_cache_integrity_one_node(node n, bool print) {
  if (n == nullptr) return;
  if (print) std::cout << " ... checking cache integrity for node " << n->key << std::endl;

  // debug check : redo the linear calculation
  auto &ca = n->cache;
  for (int b = 0; b < n_blocks; ++b) {
    auto check = check_one_block_table_linear(n, b, false);
    if (ca.block_table[b] != check) {
      std::cout << " Inconsistent block table for block " << b << " : cache =  " << ca.block_table[b] << " while it should be  " << check
                << std::endl;
      check_one_block_table_linear(n, b, true);
      TRIQS_RUNTIME_ERROR << " FATAL ";
    }
  }
}
