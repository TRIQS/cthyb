#define TRIQS_EXCEPTION_SHOW_CPP_TRACE
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#include <triqs/hilbert_space/state.hpp>
#include <triqs/hilbert_space/space_partition.hpp>
#include <triqs/hilbert_space/imperative_operator.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include "../c++/atomic_problem.hpp"
#include <iomanip>

using triqs::utility::many_body_operator;
using triqs::utility::c;
using triqs::utility::c_dag;
using triqs::hilbert_space::fundamental_operator_set;
using triqs::hilbert_space::hilbert_space;
using triqs::hilbert_space::state;
using triqs::hilbert_space::imperative_operator;
using cthyb::atomic_problem;
using triqs::arrays::matrix;
using triqs::h5::file;

using block_matrix_t = triqs::arrays::vector<matrix<double>>;
using state_t = state<hilbert_space, double, false>;
using imperative_op_t = imperative_operator<hilbert_space, double> ;

// for two matrices
inline double dot_product(matrix<double> const &a, matrix<double> const & b) { 
  double r = 0;
  int dim1 = first_dim(a), dim2 = second_dim(a);
  if ((dim1 != second_dim(b)) || ( dim2 != first_dim(b))) TRIQS_RUNTIME_ERROR << "dot_product of matrices : size mismatch";
  for (int i = 0; i< dim1; ++i) 
    for (int j = 0; j< dim2; ++j) r+=a(i,j)*b(j,i);
  return r;
} 

// For two block matrices A, B, this is Tr (tA B)
inline double dot_product(block_matrix_t const& a, block_matrix_t const& b) {
 double r = 0;
 int n_blocks = a.size();
 for (int bl = 0; bl < n_blocks; ++bl) r += dot_product(a[bl], b[bl]);
 return r;
}

// the calculation itself with its data
struct calc {

 many_body_operator<double> h_loc;
 fundamental_operator_set fops;
 atomic_problem h_diag;
 int n_blocks;
 hilbert_space full_hs;
 block_matrix_t density_matrix;

 // Read info from file
 calc() {
  h5_read(file("density_matrix.h5", 'r'), "density_matrix", density_matrix);
  h5_read(file("h_loc.h5", 'r'), "h_loc", h_loc, fops);
  h_diag = atomic_problem{h_loc, fops};
  full_hs = hilbert_space(fops);
  n_blocks = h_diag.n_subspaces();

  std::cout << "Content of fops" << std::endl;
  for (auto const& x : fops) std::cout << x.index[0] << std::endl;

  std::cout << "Block energies " << std::endl;
  std::cout << std::setprecision(15);
  for (int i = 0; i < h_diag.get_eigensystems().size(); ++i)
   std::cout << i << '\t' << h_diag.get_block_dim(i) << '\t' << h_diag.get_eigensystems()[i].eigenvalues[0] << std::endl;
 }

 // For an operator op, make the block matrix of its elements between the blocks.
 block_matrix_t matrix_elements(many_body_operator<double> const& op) const {
  imperative_op_t imp(op, fops);
  block_matrix_t R(n_blocks);
  for (int bl = 0; bl < n_blocks; ++bl) {
   R[bl] = h_diag.make_op_matrix(imp, bl, bl).first;
   // this is in the Fock basis. Now we turn it into the eigenbasis of H.
   auto const& U = h_diag.get_eigensystems()[bl].unitary_matrix;
   R[bl] = U.transpose() * R[bl] * U; // OR THE OPPOSITE ?
  }
  return R;
 }

 // Returns the static average of the operator
 double average(many_body_operator<double> const& op) const { return dot_product(matrix_elements(op), density_matrix); }

 // From an operator, return the state made by acting with it on the vacuum
 state_t act_on_vacuum(many_body_operator<double> const& op) const {
  imperative_op_t imp(op, fops);
  // create the vacuum state in st and act on it 
  state_t st(full_hs);
  st(st.get_hilbert().get_state_index(0)) = 1.0;
  return imp(st);
 }

 // Given a state in the FULL hilbertspace, compute the average of the projector onto this state
 double average_projector(state_t const& st) const {
  // loop over all eigenvectors in all sub hilbert spaces
  // and compute the overlap with the state st
  std::vector<double> overlap(h_diag.space().size());
  int index = 0;
  for (auto& es : h_diag.get_eigensystems()) {
   for (auto& ev : es.eigenstates) {
    // promote ev to a state of the full hilbert_space
    state_t eigenstate(full_hs);
    for (long i = 0; i < ev.size(); ++i) {
     long j = full_hs.get_state_index(ev.get_hilbert().get_fock_state(i));
     eigenstate(j) = ev(i);
    }
    overlap[index++] = dot_product(st, eigenstate);
   }
  }
  if (index != h_diag.space().size()) TRIQS_RUNTIME_ERROR << "Internal error";
  // now we close the sum. The density_matrix is Block matrix, need to loop over blocks
  // while overlap are in the full hilbert space
  double r = 0;
  for (int bl = 0; bl < n_blocks; ++bl) {
   int dim = h_diag.get_block_dim(bl);
   int sh = h_diag.get_first_eigstate_of_block(bl);
   for (int a = 0; a < dim; ++a)
    for (int b = 0; b < dim; ++b)   r += overlap[sh + a] * density_matrix[bl](a, b) * overlap[sh + b];
  }
  return r;
 }
 
};
