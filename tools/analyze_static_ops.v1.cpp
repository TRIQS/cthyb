#include <triqs/hilbert_space/state.hpp>
#include <triqs/hilbert_space/space_partition.hpp>
#include <triqs/hilbert_space/imperative_operator.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include "../c++/atomic_problem.hpp"

using triqs::utility::many_body_operator;
using triqs::utility::c;
using triqs::utility::c_dag;
using triqs::hilbert_space::fundamental_operator_set;
using triqs::hilbert_space::hilbert_space;
using triqs::hilbert_space::state;
using triqs::hilbert_space::imperative_operator;
using cthyb::atomic_problem;

int main() {

  // Read info from file
  many_body_operator<double> h_loc;
  fundamental_operator_set fops;
  auto f = triqs::h5::file("h_loc.h5", 'r');
  h5_read(f, "h_loc", h_loc, fops);

  std::cout << "Content of fops" << std::endl;
  for (auto const & x: fops) std::cout << x.index[0] << std::endl;

  // Part I: construct our favourite state to project on
  hilbert_space full_hs(fops);
  auto rvb_op = c_dag("00-up",0) * c_dag("01-up",0) * c_dag("11-up",0) + c_dag("01-down",0);
  imperative_operator<hilbert_space, double> rvb_imp(rvb_op, fops);

  // create the vacuum state in st and act on it with rvb_imp
  state<hilbert_space, double, false> st(full_hs);
  st(st.get_hilbert().get_state_index(0)) = 1.0;
  auto rvb_state = rvb_imp(st);

  std::cout << "vacuum state: " << st << std::endl;
  std::cout << "RVB state: " << rvb_state << std::endl;

  // Part II: let's make sure rvb_state is decomposed in a single sub hilbert space
  atomic_problem h_diag{h_loc, fops};
  int subspace_number = -1;
  for(long i=0; i<rvb_state.size(); ++i) {
    if(std::abs(rvb_state(i)) > 1.e-10) {
      long j = full_hs.get_fock_state(i); // of course j=i...
      // look for the fock state in all subspaces
      for(int k=0; k<h_diag.n_subspaces(); k++) {
        if(h_diag.subspace(k).has_state(j)) {
          if (subspace_number == -1) subspace_number = k;
          else if (subspace_number != k) std::cout << "We have a problem!!! k = " << k << " and subspace_number = " << subspace_number << std::endl;
        }
      }
    }
  }

  // Part III: get the eigenvectors and promote them to
  // the full hilbert space
  // loop over all eigenvectors in all sub hilbert spaces
  for(auto & es: h_diag.get_eigensystems()) {
    for(auto & ev: es.eigenstates) {

      state<hilbert_space, double, false> eigenstate(full_hs);
      for (long i = 0; i < ev.size(); ++i) {
        long j = full_hs.get_state_index(ev.get_hilbert().get_fock_state(i));
        eigenstate(j) = ev(i);
      }

      auto d = dot_product(rvb_state, eigenstate);
      if (std::abs(d) > 1.e-10) {
        std::cout << "dot product is: " << dot_product(rvb_state, eigenstate) << std::endl;
        std::cout << "eigenstate is: " << eigenstate << std::endl;
      }
    }
  }

}
