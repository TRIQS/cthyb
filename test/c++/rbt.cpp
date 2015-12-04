#include <triqs/utility/rbt.hpp>
#include <triqs/test_tools/arrays.hpp>
#include <iostream>
#include <fstream>
#include <map>

// We need this dummy wrapper type, because rb_tree is inherited
// from its second template parameter, and C++ forbids inheritance
// from fundamental types.
struct int_ {
  int i;
  int_(int i = {}) : i(i) {}
};

TEST(CtHyb, RBT) {

 triqs::utility::rb_tree<int,int_> tree;

 auto plot = [](triqs::utility::rb_tree<int,int_> & tree) {
 tree.graphviz(std::cout);
 tree.clear_modified();
 };

 tree.insert(1,1);
 tree.clear_modified();

 tree.insert(7,7);
 plot(tree);

 tree.insert(3,3);
 tree.insert(5,5);
 tree.insert(0,0);
 tree.insert(2,2);

 tree.insert(10,0);
 tree.clear_modified();
 tree.insert(12,2);
 plot(tree);

 tree.delete_node(7);
 plot(tree);

 tree.insert(7,7);
 plot(tree);

 tree.insert(6,6);
 plot(tree);

 triqs::utility::rb_tree<int,int_> tree_copy(tree);

 tree.delete_node(10);
 plot(tree);

 tree.delete_node(1);
 tree.delete_node(3);
 tree.delete_node(5);
 //tree.delete_node(10);
 tree.delete_node(12);
 tree.delete_node(2);
 plot(tree);

 plot(tree_copy);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
