#pragma once
#include <triqs/utility/first_include.hpp>
#include <triqs/utility/exceptions.hpp>
#include <limits>

struct rbt_insert_error {};

template <typename Key, typename Value, typename Compare = std::less<Key>> class rb_tree {

 static const bool RED = true;
 static const bool BLACK = false;
 Compare compare;

 public:
 struct node_t;
 using node = node_t*;

 struct node_t : public Value {
  Key key;          // key
                    //  Value val;        // associated data
  bool color;       // color of parent link
  int N;            // subtree count
  node left, right; // links to left and right subtrees
  bool modified, soft_deleted;

  node_t(Key const& key, Value const& val, bool color, int N)
      //     : key(key), val(val), color(color), N(N), left{nullptr}, right{nullptr}, modified(true), soft_deleted(false) {}
      : Value(val),
        key(key),
        color(color),
        N(N),
        left{nullptr},
        right{nullptr},
        modified(true),
        soft_deleted(false) {}

  node_t(node_t const&) = delete;
  node_t& operator=(node_t const&) = delete;
  template <typename... T> void reset(Key const& k, T&&... x) {
   key = k;
   left = nullptr; 
   right = nullptr;
   Value::reset(std::forward<T>(x)...);
  }
 }; 

 private:
 node root; // root of the BST

 template <typename Fnt> void apply_recursive(Fnt const& f, node n) const {
  if (n->left) apply_recursive(f, n->left);
  f(n);
  if (n->right) apply_recursive(f, n->right);
 }

 template <typename Fnt> friend void foreach(rb_tree const & tr, Fnt const& f) { if (tr.root) tr.apply_recursive(f, tr.root); }
 template <typename Fnt> friend void foreach(rb_tree const & tr, node n, Fnt const& f) { if (n) tr.apply_recursive(f, n); }

 template <typename Fnt> void apply_recursive_reverse(Fnt const& f, node n) const {
  if (n->right) apply_recursive_reverse(f, n->right);
  f(n);
  if (n->left) apply_recursive_reverse(f, n->left);
 }

 template <typename Fnt> friend void foreach_reverse(rb_tree const & tr, Fnt const& f) { if (tr.root) tr.apply_recursive_reverse(f, tr.root); }
 template <typename Fnt> friend void foreach_reverse(rb_tree const & tr, node n, Fnt const& f) { if (n) tr.apply_recursive_reverse(f, n); }

 template <typename Fnt> void apply_recursive_subtree_first(Fnt const& f, node n) const {
  if (n->left) apply_recursive_subtree_first(f, n->left);
  if (n->right) apply_recursive_subtree_first(f, n->right);
  f(n);
 }

 template <typename Fnt> friend void foreach_subtree_first(rb_tree const & tr, node n, Fnt const& f) { if (n) tr.apply_recursive_subtree_first(f, n); }
 template <typename Fnt> friend void foreach_subtree_first(rb_tree const & tr,  Fnt const& f) { foreach_subtree_first(tr,tr.root,f);}

 // is node x red; false if x is nullptr ?
 bool is_red(node x) {
  if (x == nullptr) return false;
  return (x->color == RED);
 }

 // number of node in subtree rooted at x; 0 if x is nullptr
 int size(node x) const {
  if (x == nullptr) return 0;
  return x->N;
 }

 void rec_free(node n) {
  if (n == nullptr) return;
  rec_free(n->left);
  rec_free(n->right);
  delete n;
 }

 public:
 rb_tree() : root(nullptr) {}
 ~rb_tree() { rec_free(root); }

 rb_tree(rb_tree const &) = delete; // to do 

 int size() const { return size(root); }
 bool empty() const { return root == nullptr; }
 node const &get_root() const { return root; }
 node & get_root() { return root; }

 Compare const & comparator() const { return compare;}

 void print(std::ostream & out) const {
 apply_recursive([&out](node n) {out << n->key << std::endl; }, root);
 }

 void graphviz(std::ostream&& out) const { graphviz(out); }

 void graphviz(std::ostream& out) const {
  out << "digraph G{ graph [ordering=\"out\"];" << std::endl;
  if (root) out << double(root->key) << "[color = " << (root->modified ? "red" : "black") << "];" << std::endl;
  auto f = [&out](node n) {
   if (!n) return;
   if (n->left)
    out << double(n->left->key) << "[color = " << (n->left->modified ? "red" : "black") << "];\n" <<  double(n->key) << " -> " << double(n->left->key)
        << (n->left->color == RED ? "[color = red];" : ";") << std::endl;
   if (n->right)
    out <<  double(n->right->key) << "[color = " << (n->right->modified ? "red" : "black") << "];\n" << double(n->key) << " -> " << double(n->right->key)
        << (n->right->color == RED ? "[color = red];" : ";") << std::endl;
  };
  foreach(*this,f);
  out << "}" << std::endl;
 }

 void check_no_node_modified() const {
  foreach_subtree_first(*this, [&](node y) {
   if (y && y->modified) std::cout << "node modified " << y->key << std::endl;
  });
 }
 void check_no_node_soft_deleted() const {
  foreach_subtree_first(*this, [&](node y) {
   if (y && y->modified) std::cout << "node soft deleted  " << y->key << std::endl;
  });
 }

 int clear_modified() { 
  int r = clear_modified_impl(root); 
#ifdef TRIQS_RBT_CHECKS
  check_no_node_modified();
  check_no_node_soft_deleted();
#endif
  return r;
 }

 private:
 int clear_modified_impl(node n) {
  int r = 0;
  if (n && n->modified) {
   n->modified = false;
   r++;
   if (n->soft_deleted) TRIQS_RUNTIME_ERROR << " node " << n->key << " is soft deleted";
   n->soft_deleted = false;
   r += clear_modified_impl(n->left);
   r += clear_modified_impl(n->right);
  }
  return r;
 }

 public:

 /*************************************************************************
  *  Apply a function from root to key (or null)
  *************************************************************************/

 template <typename Fnt> void get_and_apply(Key const& key, Fnt f) const {
  get_and_apply (root,key,f);
 }
 private:
 template <typename Fnt> node get_and_apply(node x, Key const& key, Fnt const& f) const {
  while (x != nullptr) {
   f(x);
   if (compare(key, x->key))
    x = x->left;
   else if (compare(x->key, key))
    x = x->right;
   else
    return x;
  }
  return nullptr;
 }

 public:

 void set_modified_from_root_to(Key const & key) { 
  get_and_apply(key,[](node y){ y->modified=true;});
 }
 
 /*************************************************************************
  *  Standard BST search
  *************************************************************************/

 // value associated with the given key; nullptr if no such key
 node get(Key const& key) const { return get(root, key); }

  // is there a key-value pair with the given key?
 bool contains(Key const& key) const { return (get(key) != nullptr); }

 // is there a key-value pair with the given key in the subtree rooted at x?
 bool contains(node x, Key const& key) const { return (get(x, key) != nullptr); }

 // value associated with the given key in subtree rooted at x; nullptr if no such key
 private:
 node get(node x, Key const& key) const {
  while (x != nullptr) {
   if (compare(key,x->key))
    x = x->left;
   else if (compare(x->key,key))
    x = x->right;
   else
    return x;
  }
  return nullptr;
 }
 /*************************************************************************
   * find_if, traversing in an order way (traversal order is fixed).
   *************************************************************************/

 // value associated with the given key; nullptr if no such key
 template <typename Fnt> friend node find_if(rb_tree const& tr, Fnt f) { return tr.find_if_impl(tr.root, f); }

 // value associated with the given key in subtree rooted at x; nullptr if no such key
 private:
 template <typename Fnt> node find_if_impl(node x, Fnt& f) const {
  if (x == nullptr) return nullptr;
  auto r = find_if_impl(x->left, f);
  if (r) return r;
  if (f(x)) return x;
  return find_if_impl(x->right, f);
 };

/*************************************************************************
  *  Red-black insertion
  *************************************************************************/
 public:

 // insert the key-value pair; overwrite the old value with the new value
 // if the key is already present
 void insert(Key const& key, Value const& val) {
  root = insert(root, key, val);
  root->color = BLACK;
  // check();
 }

 private:
 // insert the key-value pair in the subtree rooted at h
 node insert(node h, Key const& key, Value const& val) {
  if (h == nullptr) return new node_t(key, val, true, 1);

  if (compare(key,h->key))
   h->left = insert(h->left, key, val);
  else if (compare(h->key,key))
   h->right = insert(h->right, key, val);
  else
   throw rbt_insert_error{};

  // fix-up any right-leaning links
  if (is_red(h->right) && !is_red(h->left)) h = rotateLeft(h);
  if (is_red(h->left) && is_red(h->left->left)) h = rotateRight(h);
  if (is_red(h->left) && is_red(h->right)) flipColors(h);
  h->N = size(h->left) + size(h->right) + 1;

  h->modified = true;
  return h;
 }
 
 /*************************************************************************
  *  Red-black deletion
  *************************************************************************/
 private:
 // delete the key-value pair with the minimum key rooted at h
 node deleteMin(node h) {
  if (h->left == nullptr) {
   // std::cout << " dealloc " << h.p << std::endl;
   delete h;
   return nullptr;
  }
  if (!is_red(h->left) && !is_red(h->left->left)) h = moveRedLeft(h);
  h->left = deleteMin(h->left);
  return balance(h);
 }

 // delete the key-value pair with the maximum key rooted at h
 node deleteMax(node h) {
  if (is_red(h->left)) h = rotateRight(h);
  if (h->right == nullptr) {
   // std::cout << " deleting " << h->key << std::endl;
   delete h;
   return nullptr;
  }
  if (!is_red(h->right) && !is_red(h->right->left)) h = moveRedRight(h);
  h->right = deleteMax(h->right);
  return balance(h);
 }

 public:
 // delete the key-value pair with the minimum key
 void deleteMin() {
  if (empty()) TRIQS_RUNTIME_ERROR << "BST underflow";
  // if both children of root are black, set root to red
  if (!is_red(root->left) && !is_red(root->right)) root->color = RED;
  root = deleteMin(root);
  if (!empty()) root->color = BLACK;
  check();
 }

 // delete the key-value pair with the maximum key
 void deleteMax() {
  if (empty()) TRIQS_RUNTIME_ERROR << "BST underflow";
  // if both children of root are black, set root to red
  if (!is_red(root->left) && !is_red(root->right)) root->color = RED;
  root = deleteMax(root);
  if (!empty()) root->color = BLACK;
  check();
 }

 public:
 // delete the key-value pair with the given key
 void delete_node(Key const& key) {
  if (!contains(key)) TRIQS_RUNTIME_ERROR << "symbol table does not contain " << key;
  // if both children of root are black, set root to red
  if (!is_red(root->left) && !is_red(root->right)) root->color = RED;
  root = delete_node(root, key);
  if (!empty()) root->color = BLACK;
  check();
 }

 private:
 // delete the key-value pair with the given key rooted at h
 node delete_node(node h, Key const& key) {
  if (!contains(h, key)) TRIQS_RUNTIME_ERROR << " oops";

  if (compare(key,h->key)) {
   if (!is_red(h->left) && !is_red(h->left->left)) h = moveRedLeft(h);
   h->left = delete_node(h->left, key);
  } else {

   if (is_red(h->left)) h = rotateRight(h);
   if (key == h->key && (h->right == nullptr)) {
    delete h;
    return nullptr;
   }
   if (!is_red(h->right) && !is_red(h->right->left)) h = moveRedRight(h);
   if (key == h->key) {
    node x = min(h->right);
    h->key = x->key;
    h->Value::operator=(*x);
    h->modified=true; // not sure it is needed
    h->soft_deleted = false; // CRUCIAL !
    // h->val = x->val;
    // h->val = get(h->right, min(h->right).key);
    // h->key = min(h->right).key;
    h->right = deleteMin(h->right);
   } else
    h->right = delete_node(h->right, key);
  }
  return balance(h);
 }

 /*************************************************************************
  *  red-black tree helper functions
  *************************************************************************/

 // make a left-leaning link lean to the right
 node rotateRight(node h) {
  _assert((h != nullptr) && is_red(h->left));
  node x = h->left;
  h->left = x->right;
  x->right = h;
  x->color = x->right->color;
  x->right->color = RED;
  x->N = h->N;
  h->N = size(h->left) + size(h->right) + 1;
  h->modified = true;
  x->modified = true;
  return x;
 }

 // make a right-leaning link lean to the left
 node rotateLeft(node h) {
  _assert((h != nullptr) && is_red(h->right));
  node x = h->right;
  h->right = x->left;
  x->left = h;
  x->color = x->left->color;
  x->left->color = RED;
  x->N = h->N;
  h->N = size(h->left) + size(h->right) + 1;
  h->modified = true;
  x->modified = true;
  return x;
 }

 // flip the colors of a node and its two children
 void flipColors(node h) {
  // h must have opposite color of its two children
  _assert((h != nullptr) && (h->left != nullptr) && (h->right != nullptr));
  _assert((!is_red(h) && is_red(h->left) && is_red(h->right)) || ((is_red(h) && !is_red(h->left) && !is_red(h->right))));
  h->color = !h->color;
  h->left->color = !h->left->color;
  h->right->color = !h->right->color;
 }

 // Assuming that h is red and both h->left and h->left->left
 // are black, make h->left or one of its children red.
 node moveRedLeft(node h) {
  _assert((h != nullptr));
  _assert(is_red(h) && !is_red(h->left) && !is_red(h->left->left));

  flipColors(h);
  if (is_red(h->right->left)) {
   h->right = rotateRight(h->right);
   h = rotateLeft(h);
  }
  return h;
 }

 // Assuming that h is red and both h->right and h->right->left
 // are black, make h->right or one of its children red.
 node moveRedRight(node h) {
  _assert((h != nullptr));
  _assert(is_red(h) && !is_red(h->right) && !is_red(h->right->left));
  flipColors(h);
  if (is_red(h->left->left)) {
   h = rotateRight(h);
  }
  return h;
 }

 // restore red-black tree invariant
 node balance(node h) {
  _assert((h != nullptr));

  if (is_red(h->right)) h = rotateLeft(h);
  if (is_red(h->left) && is_red(h->left->left)) h = rotateRight(h);
  if (is_red(h->left) && is_red(h->right)) flipColors(h);

  h->N = size(h->left) + size(h->right) + 1;
  h->modified = true;
  return h;
 }

 /*************************************************************************
  *  Utility functions
  *************************************************************************/

 public:
 // height of tree (1-node tree has height 0)
 int height() const { return height(root); }

 private:
 int height(node x) const {
  if (x == nullptr) return -1;
  return 1 + std::max(height(x->left), height(x->right));
 }

 /*************************************************************************
  *  Ordered symbol table methods.
  *************************************************************************/
 public:
 // the smallest key; nullptr if no such key
 Key min_key() const {
  //if (empty()) return nullptr;
  return min(root)->key;
 }

 Key min_key(node x) const { return min(x)->key; }

 // the smallest key in subtree rooted at x; nullptr if no such key
 node min(node x) const {
  _assert(x != nullptr);
  if (x->left == nullptr)
   return x;
  else
   return min(x->left);
 }

 // the largest key; nullptr if no such key
 public:
 Key max_key() const {
  //if (empty()) return nullptr;
  return max(root)->key;
 }

 Key max_key(node x) const { return max(x)->key; }

 // the largest key in the subtree rooted at x; nullptr if no such key
 node max(node x) const {
  _assert(x != nullptr);
  if (x->right == nullptr)
   return x;
  else
   return max(x->right);
 }

 // the largest key less than or equal to the given key
 public:
 Key floor(Key const& key) const {
  node x = floor(root, key);
  if (x == nullptr)
   return nullptr;
  else
   return x->key;
 }

 // the largest key in the subtree rooted at x less than or equal to the given key
 private:
 node floor(node x, Key const& key) const {
  if (x == nullptr) return nullptr;
  if (key == x->key) return x;
  if (compare(key,x->key)) return floor(x->left, key);
  node t = floor(x->right, key);
  if (t != nullptr)
   return t;
  else
   return x;
 }

 // the smallest key greater than or equal to the given key
 public:
 Key ceiling(Key const& key) const {
  node x = ceiling(root, key);
  if (x == nullptr)
   return nullptr;
  else
   return x->key;
 }

 // the smallest key in the subtree rooted at x greater than or equal to the given key
 private:
 node ceiling(node x, Key const& key) const {
  if (x == nullptr) return nullptr;
  if (key == x->key) return x;
  if (compare(x->key,key)) return ceiling(x->right, key);
  node t = ceiling(x->left, key);
  if (t != nullptr)
   return t;
  else
   return x;
 }

 // the key of rank k
 public:
 Key select(int k) const {
  if (k < 0 || k >= size()) TRIQS_RUNTIME_ERROR << " unknow key"; // return nullptr;
  node x = select(root, k);
  return x->key;
 }

 // the key of rank k in the subtree rooted at x
 private:
 node select(node x, int k) const {
  _assert(x != nullptr);
  _assert(k >= 0 && k < size(x));
  int t = size(x->left);
  if (t > k)
   return select(x->left, k);
  else if (t < k)
   return select(x->right, k - t - 1);
  else
   return x;
 }

 // number of keys less than key
 public:
 int rank(Key const& key) const { return rank(key, root); }

 // number of keys less than key in the subtree rooted at x
 private:
 int rank(Key const& key, node x) const {
  if (x == nullptr) return 0;
  if (compare(key,x->key))
   return rank(key, x->left);
  else if (compare(x->key,key))
   return 1 + size(x->left) + rank(key, x->right);
  else
   return size(x->left);
 }

 /*************************************************************************
   *  DEBUG CODE : Check integrity of red-black BST data structure
   *************************************************************************/
 private:
 void _assert(bool c) const {
  if (!c) TRIQS_RUNTIME_ERROR << "Error";
 }

 bool check() {

  return true;

  /* if (!isBST()) TRIQS_RUNTIME_ERROR << "Not in symmetric order";
   if (!isSizeConsistent()) TRIQS_RUNTIME_ERROR << "Subtree counts not consistent";
   if (!isRankConsistent()) TRIQS_RUNTIME_ERROR << "Ranks not consistent";
   if (!is23()) TRIQS_RUNTIME_ERROR << "Not a 2-3 tree";
   if (!isBalanced()) TRIQS_RUNTIME_ERROR << "Not balanced";
   return isBST() && isSizeConsistent() && isRankConsistent() && is23() && isBalanced();
 */
 }

 // does this binary tree satisfy symmetric order?
 // Note: this test also ensures that data structure is a binary tree since order is strict
 bool isBST() { return isBST(root, std::numeric_limits<int>::min(), std::numeric_limits<int>::max()); }
 // bool isBST() { return isBST(root, std::numeric_limits<int>::min(), std::numeric_limits<int>::max()); }

 // is the tree rooted at x a BST with all keys strictly between min and max
 // (if min or max is nullptr, treat as empty constraint)
 // Credit: Bob Dondero's elegant solution
 bool isBST(node x, Key min, Key max) {
  if (x == nullptr) return true;
  // if (min != nullptr && x->key <= min) return false;
  // if (max != nullptr && x->key >= max) return false;
  return isBST(x->left, min, x->key) && isBST(x->right, x->key, max);
 }

 // are the size fields correct?
 bool isSizeConsistent() { return isSizeConsistent(root); }

 bool isSizeConsistent(node x) {
  if (x == nullptr) return true;
  if (x->N != size(x->left) + size(x->right) + 1) return false;
  return isSizeConsistent(x->left) && isSizeConsistent(x->right);
 }

 // check that ranks are consistent
 bool isRankConsistent() {
  for (int i = 0; i < size(); i++)
   if (i != rank(select(i))) return false;
  // for (Key const& key : keys())
  // if (key != select(rank(key))) return false;
  return true;
 }

 // Does the tree have no red right links, and at most one (left)
 // red links in a row on any path?
 bool is23() { return is23(root); }
 bool is23(node x) {
  if (x == nullptr) return true;
  if (is_red(x->right)) return false;
  if (x != root && is_red(x) && is_red(x->left)) return false;
  return is23(x->left) && is23(x->right);
 }

 // do all paths from root to leaf have same number of black edges?
 bool isBalanced() {
  int black = 0; // number of black links on path from root to min
  node x = root;
  while (x != nullptr) {
   if (!is_red(x)) black++;
   x = x->left;
  }
  return isBalanced(root, black);
 }

 // does every path from the root to a leaf have the given number of black links?
 bool isBalanced(node x, int black) {
  if (x == nullptr) return black == 0;
  if (!is_red(x)) black--;
  return isBalanced(x->left, black) && isBalanced(x->right, black);
 }
};

