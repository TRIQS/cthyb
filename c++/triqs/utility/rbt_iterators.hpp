/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014 by O. Parcollet, M. Ferrero, P. Seth
 * Adapted from Algorithms (fourth edition) by R. Sedgewick
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
#include <stack>

namespace triqs {
  namespace utility {

    template <typename RBT> using get_node_t = std14::conditional_t<std::is_const<RBT>::value, typename RBT::node const, typename RBT::node>;

    // flatten the tree in ascending order
    template <typename RBT> std::vector<get_node_t<RBT>> flatten2(RBT &tree) {
      using node_t = get_node_t<RBT>;
      std::vector<node_t> R;
      R.reserve(tree.size());
      std::stack<node_t> stack;
      node_t n = tree.get_root();
      while (1) {
        while (n != nullptr) {
          stack.push(n);
          n = n->left;
        }
        if (stack.size() == 0) break;
        n = stack.top();
        stack.pop();
        R.push_back(n);
        n = n->right;
      }
      return R;
    }

    template <typename RBT> std::vector<get_node_t<RBT>> flatten(RBT &tree) {
      using node_t = get_node_t<RBT>;
      std::vector<node_t> R;
      R.reserve(tree.size());
      foreach (tree, [&R](node_t n) { R.push_back(n); })
        ;
      return R;
    }

    /// implementation of iterators
    namespace detail {

      // forward iterator
      template <typename RBT, typename Node> class rbt_iterator : public std::iterator<std::forward_iterator_tag, Node> {

        RBT *tree    = nullptr;
        Node n       = nullptr;
        Node current = nullptr;
        std::stack<Node> stack;

        public:
        rbt_iterator() = default;
        rbt_iterator(RBT *tree, bool at_end) : tree(tree), n(at_end ? nullptr : tree->get_root()) { operator++(); }

        rbt_iterator &operator++() {
          while (n != nullptr) {
            stack.push(n);
            n = n->left;
          }
          if (stack.size() != 0) {
            n = stack.top();
            stack.pop();
            current = n;
            n       = n->right;
          } else
            current = nullptr;
          return *this;
        }

        Node &operator*() { return current; }
        Node &operator->() { return current; }

        rbt_iterator operator++(int) {
          auto c = *this;
          operator++();
          return c;
        }

        bool operator==(rbt_iterator const &other) const { return (other.current == current); }
        bool operator!=(rbt_iterator const &other) const { return (!operator==(other)); }
      };
    }
  }
}
