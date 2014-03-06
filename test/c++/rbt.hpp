#pragma once
#include <triqs/utility/first_include.hpp>
#include <triqs/utility/exceptions.hpp>
#include <limits>

template <typename Key, typename Value> class rb_tree {

 static const bool RED = true;
 static const bool BLACK = false;

 struct node_t;
 using node = node_t *; 

 struct node_t  { 
  Key key;          // key
  Value val;        // associated data
  bool color;       // color of parent link
  int N;            // subtree count
  node left, right; // links to left and right subtrees
  bool modified;    

  node_t(Key const& key, Value const& val, bool color, int N)
     : key(key), val(val), color(color), N(N), left{nullptr}, right{nullptr}, modified(true) {}
 }; // node_t

 node root; // root of the BST

 template <typename Fnt> void apply_recursive(Fnt const& f, node n) const {
  if (n->left) apply_recursive(f, n->left);
  f(n);
  if (n->right) apply_recursive(f, n->right);
 }
 template <typename Fnt> void apply_recursive(Fnt const& f) const { apply_recursive(f, root); }

 // is node x red; false if x is nullptr ?
 bool is_red(node x) {
  if (x == nullptr) return false;
  return (x->color == RED);
 }

 // number of node in subtree rooted at x; 0 if x is nullptr
 int size(node x) {
  if (x == nullptr) return 0;
  return x->N;
 }

 void rec_free(node n) {
  if (n==nullptr) return;
  rec_free(n->left);
  rec_free(n->right);
  delete n;
 }

 public:
 rb_tree() : root(nullptr) {}

 ~rb_tree() { rec_free(root); }

 int size() { return size(root); }
 bool empty() { return root == nullptr; }

 //void print(std::ostream & out) const {
 // apply_recursive([](node n) {out << n->key << std::endl; });
 //}

 void graphviz(std::ostream && out) const { graphviz(out);}

 void graphviz(std::ostream & out) const {
  out << "digraph G{ graph [ordering=\"out\"];" << std::endl;
  if (root) out << root->key << "[color = " << (root->modified ? "red" : "black") << "];"<< std::endl;
  auto f = [&out](node n) {
   if (!n) return;
   if (n->left)
    out << n->left->key << "[color = " << (n->left->modified ? "red" : "black") << "];\n" << n->key << " -> " << n->left->key
        << (n->left->color == RED ? "[fillcolor = red];" : ";") << std::endl;
   if (n->right) 
       out << n->right->key << "[color = " << (n->right->modified ? "red" : "black") << "];\n" << n->key << " -> " << n->right->key
        << (n->right->color == RED ? "[fillcolor = red];" : ";") << std::endl;
  };
  apply_recursive(f);
  out << "}" << std::endl;
 }


 void clear_modified_impl (node n) { 
  if (n && n->modified) {
   n->modified = false; 
   clear_modified_impl (n->left);
   clear_modified_impl (n->right);
  }
 }
 
 void clear_modified () { clear_modified_impl(root);} 
 
 /*************************************************************************
  *  Standard BST search
  *************************************************************************/

 // value associated with the given key; nullptr if no such key
 node get(Key const& key) { return get(root, key); }

 // value associated with the given key in subtree rooted at x; nullptr if no such key
 private:
 node get(node x, Key const& key) {
  while (x != nullptr) {
   if (key < x->key)
    x = x->left;
   else if (key > x->key)
    x = x->right;
   else
    return x;
  }
  return nullptr;
 }

 public:
 // is there a key-value pair with the given key?
 bool contains(Key const& key) { return (get(key) != nullptr); }

 // is there a key-value pair with the given key in the subtree rooted at x?
 bool contains(node x, Key const& key) { return (get(x, key) != nullptr); }

 /*************************************************************************
  *  Red-black insertion
  *************************************************************************/

 // insert the key-value pair; overwrite the old value with the new value
 // if the key is already present
 void insert(Key const& key, Value const& val) {
  root = insert(root, key, val);
  root->color = BLACK;
  check();
 }

 private:
 // insert the key-value pair in the subtree rooted at h
 node insert(node h, Key const& key, Value const& val) {
  if (h == nullptr) return new node_t(key, val, true, 1);

  if (key < h->key)
   h->left = insert(h->left, key, val);
  else if (key > h->key)
   h->right = insert(h->right, key, val);
  else
   h->val = val; // ERROR TREATMENT

  // RECOMPUTE 

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
   //std::cout << " dealloc " << h.p << std::endl;
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
   //std::cout << " dealloc " << h.p << std::endl;
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
  // if (!contains(h, key)) TRIQS_RUNTIME_ERROR << " oops";

  if (key < h->key) {
   if (!is_red(h->left) && !is_red(h->left->left)) h = moveRedLeft(h);
   h->left = delete_node(h->left, key);
  } else {

   if (is_red(h->left)) h = rotateRight(h);
   if (key == h->key && (h->right == nullptr)) {
    //std::cout << " dealloc xx : " << h.p << std::endl;
    delete h;
    return nullptr;
   }
   if (!is_red(h->right) && !is_red(h->right->left)) h = moveRedRight(h);
   if (key == h->key) {
    node x = min(h->right);
    h->key = x->key;
    h->val = x->val;
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
 int height() { return height(root); }

 private:
 int height(node x) {
  if (x == nullptr) return -1;
  return 1 + std::max(height(x->left), height(x->right));
 }

 /*************************************************************************
  *  Ordered symbol table methods.
  *************************************************************************/
 public:
 // the smallest key; nullptr if no such key
 Key min() {
  if (empty()) return nullptr;
  return min(root).key;
 }

 // the smallest key in subtree rooted at x; nullptr if no such key
 private:
 node min(node x) {
  _assert(x != nullptr);
  if (x->left == nullptr)
   return x;
  else
   return min(x->left);
 }

 // the largest key; nullptr if no such key
 public:
 Key max() {
  if (empty()) return nullptr;
  return max(root).key;
 }

 // the largest key in the subtree rooted at x; nullptr if no such key
 private:
 node max(node x) {
  _assert(x != nullptr);
  if (x->right == nullptr)
   return x;
  else
   return max(x->right);
 }

 // the largest key less than or equal to the given key
 public:
 Key floor(Key const& key) {
  node x = floor(root, key);
  if (x == nullptr)
   return nullptr;
  else
   return x->key;
 }

 // the largest key in the subtree rooted at x less than or equal to the given key
 private:
 node floor(node x, Key const& key) {
  if (x == nullptr) return nullptr;
  if (key == x->key) return x;
  if (key < x->key) return floor(x->left, key);
  node t = floor(x->right, key);
  if (t != nullptr)
   return t;
  else
   return x;
 }

 // the smallest key greater than or equal to the given key
 public:
 Key ceiling(Key const& key) {
  node x = ceiling(root, key);
  if (x == nullptr)
   return nullptr;
  else
   return x->key;
 }

 // the smallest key in the subtree rooted at x greater than or equal to the given key
 private:
 node ceiling(node x, Key const& key) {
  if (x == nullptr) return nullptr;
  if (key == x->key) return x;
  if (key > x->key) return ceiling(x->right, key);
  node t = ceiling(x->left, key);
  if (t != nullptr)
   return t;
  else
   return x;
 }

 // the key of rank k
 public:
 Key select(int k) {
  if (k < 0 || k >= size()) TRIQS_RUNTIME_ERROR << " unknow key"; // return nullptr;
  node x = select(root, k);
  return x->key;
 }

 // the key of rank k in the subtree rooted at x
 private:
 node select(node x, int k) {
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
 int rank(Key const& key) { return rank(key, root); }

 // number of keys less than key in the subtree rooted at x
 private:
 int rank(Key const& key, node x) {
  if (x == nullptr) return 0;
  if (key < x->key)
   return rank(key, x->left);
  else if (key > x->key)
   return 1 + size(x->left) + rank(key, x->right);
  else
   return size(x->left);
 }

 /*************************************************************************
   *  DEBUG CODE : Check integrity of red-black BST data structure
   *************************************************************************/
 private:
 void _assert(bool c) {
  if (!c) TRIQS_RUNTIME_ERROR << "Error";
 }

 bool check() {
  
  return true;

  if (!isBST()) TRIQS_RUNTIME_ERROR << "Not in symmetric order";
  if (!isSizeConsistent()) TRIQS_RUNTIME_ERROR << "Subtree counts not consistent";
  if (!isRankConsistent()) TRIQS_RUNTIME_ERROR << "Ranks not consistent";
  if (!is23()) TRIQS_RUNTIME_ERROR << "Not a 2-3 tree";
  if (!isBalanced()) TRIQS_RUNTIME_ERROR << "Not balanced";
  return isBST() && isSizeConsistent() && isRankConsistent() && is23() && isBalanced();
 }

 // does this binary tree satisfy symmetric order?
 // Note: this test also ensures that data structure is a binary tree since order is strict
 bool isBST() { return isBST(root, std::numeric_limits<int>::min(), std::numeric_limits<int>::max()); }
 //bool isBST() { return isBST(root, std::numeric_limits<int>::min(), std::numeric_limits<int>::max()); }

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

