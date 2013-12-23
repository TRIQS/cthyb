#pragma once
#include "krylov_worker.hpp"
#include "sorted_spaces.hpp"

namespace cthyb_krylov {

template <typename HamiltonianType, typename StateType> class exp_h_worker {

 using state_type = StateType;
 using scalar_type = typename state_type::scalar_t;

 krylov_worker<HamiltonianType, StateType> kw;
 sorted_spaces const& sosp;

 int small_matrix_size;

 // Temporary matrices
 matrix<scalar_type> matrix_exp;
 matrix<scalar_type> krylov_exp;
 triqs::arrays::vector<scalar_type> v2;

 public:
 exp_h_worker(HamiltonianType const& H, sorted_spaces const& sosp, double gs_energy_convergence, int small_matrix_size)
    : kw(H, gs_energy_convergence), sosp(sosp), small_matrix_size(small_matrix_size) {
  int max_subspace_dim = 0;
  for (int nsp = 0; nsp < sosp.n_subspaces(); ++nsp)
   max_subspace_dim = std::max(max_subspace_dim, sosp.subspace(nsp).dimension());
  small_matrix_size = std::min(small_matrix_size, max_subspace_dim);
  matrix_exp.resize(small_matrix_size, small_matrix_size);
  krylov_exp.resize(max_subspace_dim, max_subspace_dim);
  v2.resize(max_subspace_dim);
 }

 //---------------------------------------------------------------------------

 state_type operator()(state_type const& initial_state, double dtau) {
  auto const& space = initial_state.get_hilbert();
  int space_dim = space.dimension();

  if (space_dim > small_matrix_size) {

   scalar_type initial_state_norm = std::sqrt(dot_product(initial_state, initial_state));
   kw(initial_state / initial_state_norm);

   auto eigenvalues = kw.values();
   int krylov_dim = eigenvalues.size();

   auto all = range(0, krylov_dim);

   // TO BE OPTIMIZED !! Like below
   krylov_exp(all, all) = 0;
   for (int n = 0; n < krylov_dim; ++n) krylov_exp(n, n) = exp(-dtau * eigenvalues(n));

   krylov_exp(all, all) = kw.vectors().transpose() * krylov_exp(all, all) * kw.vectors();
   auto krylov_coeffs = initial_state_norm * krylov_exp(all, 0);

   return kw.krylov_2_fock(krylov_coeffs);

  } else {

   auto const& eigensystem = sosp.get_eigensystems()[space.get_index()];
   auto all = range(0, space_dim);
   StateType st = make_zero_state(initial_state);

   // v2(all) = unitary_matrix.transpose() * initial_state.amplitudes();
   triqs::arrays::blas::gemv(1.0, eigensystem.unitary_matrix.transpose(), initial_state.amplitudes(), 0.0, v2(all));
   for (int n = 0; n < space_dim; ++n) v2[n] *= exp(-dtau * eigensystem.eigenvalues(n));
   // st.amplitudes() =  unitary_matrix * v2(all);
   triqs::arrays::blas::gemv(1.0, eigensystem.unitary_matrix, v2(all), 0.0, st.amplitudes());

   return st;
  }
 }
 //-----------------------------------------------------

 void apply_no_emin(state_type& initial_state, double dtau) {
  auto const& space = initial_state.get_hilbert();
  int space_dim = space.dimension();

  auto const& eigensystem = sosp.get_eigensystems()[space.get_index()];

  if (space_dim == 1) return;
  auto all = range(0, space_dim);
  triqs::arrays::blas::gemv(1.0, eigensystem.unitary_matrix.transpose(), initial_state.amplitudes(), 0.0, v2(all));
  for (int n = 1; n < space_dim; ++n) v2[n] *= exp(-dtau * (eigensystem.eigenvalues(n) - eigensystem.eigenvalues(0)));
  triqs::arrays::blas::gemv(1.0, eigensystem.unitary_matrix, v2(all), 0.0, initial_state.amplitudes());
 }
};
}
