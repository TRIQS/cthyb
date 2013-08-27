#include "./fundamental_operator_set.hpp"
#include "./fock_state.hpp"
#include "./complete_hilbert_space.hpp"
#include "./partial_hilbert_space.hpp"
#include "./operator.hpp"
#include "./imperative_operator.hpp"
#include "./state.hpp"

using namespace triqs::app::impurity_solvers::ctqmc_krylov;

int main() {

  std::cout << std::endl << "Part I: the fundamental_operator_set class" << std::endl << std::endl;

  fundamental_operator_set<int,int> f1(std::vector<int>(2,4));
  for (int i=0; i<2; ++i)
    for (int j=0; j<4; ++j)
      std::cout << "(" << i << "," << j << ") --> " << f1.index_to_n(i,j) << std::endl;
  std::cout << "dim = " << f1.dimension() << std::endl;
  std::cout << "n operators = " << f1.n_operators() << std::endl;

  std::cout << std::endl;
  fundamental_operator_set<int,int> f2;
  f2 = f1;
  std::cout << "(1,1) --> " << f2.index_to_n(1,1) << std::endl;

  std::cout << std::endl;
  fundamental_operator_set<int> f3;
  for (int i=0; i<4; ++i) f3.add_operator(i);
  std::cout << "2 --> " << f3.index_to_n(2) << std::endl;
  std::cout << "dim = " << f3.dimension() << std::endl;
  std::cout << "n operators = " << f3.n_operators() << std::endl;

  std::cout << std::endl;
  fundamental_operator_set<std::string,int> f4;
  for (int i=0; i<2; ++i) f4.add_operator("up",i);
  for (int i=0; i<2; ++i) f4.add_operator("down",i);
  std::cout << "(down,0) --> " << f4.index_to_n("down",0) << std::endl;
  std::cout << "dim = " << f4.dimension() << std::endl;
  std::cout << "n operators = " << f4.n_operators() << std::endl;
  
  std::cout << std::endl << "Part II: the complete_hilbert_space class" << std::endl << std::endl;

  complete_hilbert_space h(f1);
  std::cout << "dim = " << h.dimension() << std::endl;
  std::cout << "fock state for index 120 = " << h.get_fock_state(120) << std::endl;
  std::cout << "index of fock state 120 = " << h.get_state_index(h.get_fock_state(120)) << std::endl;
  complete_hilbert_space h2;
  h2 = h;
  std::cout << "dim = " << h.dimension() << std::endl;
  std::cout << std::endl;
  
  std::cout << std::endl << "Part III: the fock_state class" << std::endl << std::endl;
  
  complete_hilbert_space h3(f3);
  
  fock_state fs1 = h3.get_fock_state(10);
  std::cout << fs1 << std::endl;
  fock_state fs2 = fs1;
  std::cout << fs2 << std::endl;
  std::cout << std::endl;

  std::cout << std::endl << "Part IV: the state operator" << std::endl << std::endl;

  fundamental_operator_set<const char *, int> fop;
  for (int i=0; i<5; ++i) fop.add_operator("up",i);

  complete_hilbert_space h_full(fop);
  state<complete_hilbert_space, true> st(h_full);
  st(0) = 3.0;
  st(3) = 5.0;
  std::cout << "state is: " << st << std::endl;

  std::cout << std::endl << "Part V: the declarative operator" << std::endl << std::endl;

  using triqs::utility::c;
  using triqs::utility::c_dag;
  
  auto H = 3 * c_dag("up",1) * c("up",1) + 2 * c_dag("up",2) * c("up",2) + c("up",1) * c("up",2);
  std::cout << "H = " << H << std::endl;

  std::cout << std::endl << "Part VI: the imperative operator" << std::endl << std::endl;

  auto opH = imperative_operator<complete_hilbert_space>(H, fop);

  state<complete_hilbert_space, true> old_state(h_full);
  old_state(7) = 1.0;
  std::cout << "old state is: " << old_state << std::endl;

  auto new_state = opH(old_state);
  std::cout << "new state is: " << new_state << std::endl;

  std::unordered_map<const complete_hilbert_space*, const complete_hilbert_space*> mymap;
  mymap[&h_full] = &h_full;
  auto opH2 = imperative_operator<complete_hilbert_space, true>(H, fop, mymap);

  state<complete_hilbert_space, true> old_state2(h_full);
  old_state2(7) = 1.0;
  std::cout << "old state is: " << old_state2 << std::endl;

  auto new_state2 = opH2(old_state2);
  std::cout << "new state is: " << new_state2 << std::endl;

  auto copy_op = opH2;
  std::cout << "new state is: " << copy_op(old_state2) << std::endl;

  std::cout << "size of hilbert: " << sizeof(h_full) << std::endl;
  std::cout << "size of op: " << sizeof(opH) << std::endl;
  std::cout << "size of op: " << sizeof(opH2) << std::endl;
  std::cout << "size of op: " << sizeof(copy_op) << std::endl;


  std::cout << std::endl << "Part VII: partial hilbert space" << std::endl << std::endl;
  
  auto Cdag = c_dag("up",2);

  fundamental_operator_set<const char *, int> fop2;
  for (int i=0; i<5; ++i) fop2.add_operator("up",i);
  
  complete_hilbert_space h4(f4);
  
  partial_hilbert_space phs0;
  phs0.add_basis_fock(h4.get_fock_state(0)); // 000
  phs0.add_basis_fock(h4.get_fock_state(1)); // 001
  phs0.add_basis_fock(h4.get_fock_state(2)); // 010
  phs0.add_basis_fock(h4.get_fock_state(3)); // 011

  partial_hilbert_space phs1;
  phs1.add_basis_fock(h4.get_fock_state(4)); // 100
  phs1.add_basis_fock(h4.get_fock_state(5)); // 101
  phs1.add_basis_fock(h4.get_fock_state(6)); // 110
  phs1.add_basis_fock(h4.get_fock_state(7)); // 111

  std::unordered_map<const partial_hilbert_space*, const partial_hilbert_space*> Cdagmap;
  Cdagmap[&phs0] = &phs1;
  auto opCdag = imperative_operator<partial_hilbert_space, true>(Cdag, fop2, Cdagmap);

  state<partial_hilbert_space, false> start(phs0);

  start(0) = 1.0;
  start(1) = 2.0;
  start(2) = 3.0;
  start(3) = 4.0;

  std::cout << "new state is: " << opCdag(start) << std::endl;

  return 0;

}
