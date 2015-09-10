#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/state.hpp>
#include "../c++/sorted_spaces.hpp"

using triqs::utility::many_body_operator;
using triqs::hilbert_space::fundamental_operator_set;
using triqs::hilbert_space::hilbert_space;
using triqs::hilbert_space::state;
using cthyb::sorted_spaces;
using triqs::arrays::matrix;

using block_matrix_t = std::vector<matrix<double>>;
using state_t = state<triqs::hilbert_space::hilbert_space, double, false>;

class post_process {

 many_body_operator<double> h_loc;
 fundamental_operator_set fops;
 sorted_spaces sosp;
 int n_blocks;
 triqs::hilbert_space::hilbert_space full_hs;
 block_matrix_t density_matrix;

 public:
 /// From the density_matrix and the h_loc Hamiltonian
 post_process(block_matrix_t const& density_matrix, many_body_operator<double> const& h_loc);

 /// From an operator op, make the block matrix of its elements between the blocks.
 block_matrix_t matrix_elements(many_body_operator<double> const& op) const;

 /// Returns the static average of the operator
 double average(many_body_operator<double> const& op) const;

 // From an operator, return the state made by acting with it on the vacuum
 TRIQS_WRAP_IGNORE state_t act_on_vacuum(many_body_operator<double> const& op) const;

 // Given a state in the FULL hilbertspace, compute the average of the projector onto this state
 TRIQS_WRAP_IGNORE double average_projector(state_t const& st) const;

 /// static average of the projector on the state create by op on the vaccum
 double average_projector_from_creation_op(many_body_operator<double> const& op) const {
  return average_projector(act_on_vacuum(op));
 }

 /// dot product of the 2 states created by op1 and op2 on the vacuum
 double dot_product_from_creation_op(many_body_operator<double> const& op1, many_body_operator<double> const& op2);
};
