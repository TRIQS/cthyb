#include <triqs/utility/first_include.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>

using triqs::operators::c;
using triqs::operators::c_dag;

using atom_diag= triqs::atom_diag::atom_diag<false>;
using namespace triqs::atom_diag;

int main() {

  // define operators
  auto n_up   = c_dag("up") * c("up");
  auto n_down = c_dag("down") * c("down");

  // choose the fundamental operator set
  fundamental_operator_set fops;
  fops.insert("up");
  fops.insert("down");

  // Hubbard atom in magnetic field
  {
    double U = 2.0;
    double h = 0.5;
    double d = 0.5;

    auto H = U * n_up * n_down + h * (n_up - n_down);

    // put quantum numbers in a vector
    std::vector<many_body_op_t> qn_list{n_up, n_down};

    // Divide the full Hilbert space
    atom_diag h_diag(H, fops, qn_list);
    std::cout << h_diag << std::endl;
  }

  // Half-filled Hubbard atom + spin flips + anomalous terms
  {
    double U = 2.0;
    double J = 0.3;
    double d = 0.1;

    auto H = U * (n_up - 0.5) * (n_down - 0.5);
    H += J * (c_dag("up") * c("down") + c_dag("down") * c("up"));
    H += d * (c_dag("up") * c_dag("down") - c("up") * c("down"));

    // Divide the full Hilbert space
    atom_diag h_diag(H, fops);
    std::cout << h_diag << std::endl;
  }
}
