#include "./rbt.hpp"
#include <iostream>
#include <fstream>
#include <map>
int main() { 
 rb_tree<int,int> tree;

 auto plot = [&tree]() {
 tree.graphviz(std::ofstream("g"));
 //system("dot -Tpdf g >g.pdf && open g.pdf");
 //std::cin.ignore();
 tree.clear_modified();
 };

if (1) { 

 tree.insert(1,1);
 tree.clear_modified();

 tree.insert(7,7);
 plot();

 tree.insert(3,3);
 tree.insert(5,5);
 tree.insert(0,0);
 tree.insert(2,2);

 tree.insert(10,0);
 tree.clear_modified();
 tree.insert(12,2);

 plot();
 //tree.graphviz(std::ofstream("g"));

 tree.delete_node(7);
 plot();

 tree.insert(7,7);
 plot();

 tree.insert(6,6);
 plot();

 tree.delete_node(10);
 plot();


 tree.delete_node(1);
 tree.delete_node(3);
 tree.delete_node(5);
 //tree.delete_node(10);
 tree.delete_node(12);
 tree.delete_node(2);

}
 else { 

std::map<int,int> t2;

int N =1000;
if (1) { 
 std::cout << " mon arbre"<< std::endl;
 for(int j=0; j<N; ++j)
for(int i=0; i<N; ++i)
tree.insert(i*N+j,j);

}
else {
 std::cout << "std map "<< std::endl; 
for(int j=0; j<N; ++j)
for(int i=0; i<N; ++i)
t2.insert({i*N+j,j});
}

}

}
