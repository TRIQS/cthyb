#include <triqs/utility/first_include.hpp>
#include "./sorted_spaces.hpp"
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>

using namespace cthyb;
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
  auto qn_list = std::vector<many_body_operator<double>>{n_up, n_down};

  // chose the fundamental operator set
  fundamental_operator_set fops;
  fops.insert("up");
  fops.insert("down");

  // Divide the full Hilbert space
  sorted_spaces ss(H, qn_list, fops);
  std::cout << ss << std::endl;

  return 0;
}
