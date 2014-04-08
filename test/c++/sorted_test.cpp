#include <triqs/utility/first_include.hpp>
#include "./sorted_spaces.hpp"
#include "./operator.hpp"
#include <triqs/draft/hilbert_space_tools/fundamental_operator_set.hpp>
#include "./gf_block_structure.hpp"

using namespace cthyb_matrix;
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

  // Block structure
  std::vector<block_desc_t> block_structure;
  block_structure.push_back({"up",{{"up"}}});
  block_structure.push_back({"down",{{"down"}}});

  // divide the full Hilbert space
  sorted_spaces ss(H, qn_list, fops);
  std::cout << ss << std::endl;

  return 0;
}
