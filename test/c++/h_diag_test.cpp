#include <triqs/utility/first_include.hpp>
#include "./atom_diag.hpp"
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>

using namespace cthyb;
using triqs::operators::many_body_operator_real;
using triqs::operators::c;
using triqs::operators::c_dag;

int main() {

  // define operators
  double U = 2.0;
  double h = 0.5;
  auto n_up = c_dag("up") * c("up");
  auto n_down = c_dag("down") * c("down");
  auto H = U * n_up * n_down + h*(n_up - n_down);

  // put quantum numbers in a vector
  std::vector<many_body_operator_real> qn_list{n_up, n_down};

  // chose the fundamental operator set
  fundamental_operator_set fops;
  fops.insert("up");
  fops.insert("down");

  // Divide the full Hilbert space
  atom_diag h_diag(H, fops, qn_list);
  std::cout << h_diag << std::endl;

}
