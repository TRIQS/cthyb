#include <triqs/utility/first_include.hpp>
#include "./sorted_spaces.hpp"
#include "./operator.hpp"
#include "./fundamental_operator_set.hpp"

using namespace triqs::app::impurity_solvers::ctqmc_krylov;
using triqs::utility::many_body_operator;
using triqs::utility::c;
using triqs::utility::c_dag;

int main() {

  // define operators
  double U = 2.0;
  auto n_up = c_dag("up") * c("up");
  auto n_down = c_dag("down") * c("down");
  auto H = U * n_up * n_down;

  // put quantum numbers in a vector
#ifndef TRIQS_WORKAROUND_INTEL_COMPILER_BUGS
  std::vector<many_body_operator<double,std::string>> qn_list {n_up, n_down};
#else
  std::vector<many_body_operator<double,std::string>> qn_list; qn_list.push_back(n_up); qn_list.push_back(n_down);
#endif
  //std::vector<many_body_operator<double,std::string, int>> qn_list {n_up+n_down};
  //std::vector<many_body_operator<double,std::string, int>> qn_list;

  // chose the fundamental operator set
  fundamental_operator_set<std::string> fops;
  fops.add_operator("up");
  fops.add_operator("down");

  // Block structure
  std::vector<block_desc_t<std::string>> block_structure;
  block_structure.push_back({"up",{std::make_tuple("up")}});
  block_structure.push_back({"down",{std::make_tuple("down")}});

  // divide the full Hilbert space
  sorted_spaces ss(H, qn_list, fops, block_structure);
  std::cout << ss << std::endl;

  // get a state in sub Hilbert space 0
  auto st = ss.substate(0); st(0) = 1.0;

  // get the imperative creation operator c_dag("up",1)
  auto op = ss.get_fundamental_operator(true,0,0);

  // print connection map for that operator
  for (int n=0; n<ss.n_subspaces(); ++n) {
    std::cout << n << " --> " << ss.fundamental_operator_connect(true,0,0,n) << std::endl;
  }

  // print action on state
  std::cout << st << " --> " << op(st)<< std::endl;
  return 0;

}
