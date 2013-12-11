#include "atomic_correlators_worker.hpp"
#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>
#include <algorithm>

namespace cthyb_krylov {

atomic_correlators_worker::atomic_correlators_worker(configuration& c, sorted_spaces const& sosp_, double gs_energy_convergence,
                                                     int small_matrix_size)
   : config(&c),
     sosp(sosp_),
     exp_h(sosp.get_hamiltonian(), sosp, gs_energy_convergence, small_matrix_size),
     small_matrix_size(small_matrix_size) {}

//------------------------------------------------------------------------------

atomic_correlators_worker::result_t atomic_correlators_worker::operator()() {
 auto _begin = config->oplist.crbegin();
 auto _end = config->oplist.crend();
 auto last_tau = config->beta();
 int n_blocks = sosp.n_subspaces();

//#define NO_FIRST_PASS
#ifndef NO_FIRST_PASS
 // make a first pass to compute the bound for each term.
 std::vector<double> E_min_delta_tau(n_blocks, 0);
 std::vector<int> blo(n_blocks);
 for (int u = 0; u < n_blocks; ++u) blo[u] = u;

 // do the first exp
 double dtau = (_begin == _end ? config->beta() : double(_begin->first));
 for (int n = 0; n < n_blocks; ++n) {
  E_min_delta_tau[n] += dtau * sosp.get_eigensystems()[n].eigenvalues[0]; // delta_tau * E_min_of_the_block
 }
 
 for (auto it = _begin; it != _end;) { // do nothing if no operator
  auto it1 = it;
  ++it;
  dtau = (it == _end ? config->beta() : double(it->first)) - double(it1->first);
  bool one_non_zero = false;
  for (int n = 0; n < n_blocks; ++n) {
   if (blo[n] == -1) continue; // that chain has ended
   E_min_delta_tau[n] += dtau * sosp.get_eigensystems()[blo[n]].eigenvalues[0]; // delta_tau * E_min_of_the_block
   blo[n] = sosp.fundamental_operator_connect_from_linear_index(it1->second.dagger, it1->second.linear_index, blo[n]);
   one_non_zero |= (blo[n] != -1);
   // a bit slower
   //blo[n] = sosp.fundamental_operator_connect(it1->second.dagger, it1->second.block_index, it1->second.inner_index, blo[n]);
  }
  if (!one_non_zero) return 0; // quick exit, the trace is structurally 0
 }

 // Now sort the blocks
 std::vector<std::pair<double, int>> to_sort(n_blocks);
 int n_bl = 0;// the number of blocks giving non zero
 for (int n = 0; n < n_blocks; ++n)
  if (blo[n] == n) // Must return to the SAME block, or trace is 0
   to_sort[n_bl++] = std::make_pair(E_min_delta_tau[n], n);
  std::sort(to_sort.begin(), to_sort.begin() + n_bl); // sort those vector

#endif

 result_t full_trace = 0;
 double epsilon = 1.e-15;

 // To implement : regroup all the vector of the block for dgemm computation !
#ifndef NO_FIRST_PASS
 for (int bl = 0; ((bl < n_bl) && (std::exp(- to_sort[bl].first) >= ( std::abs(full_trace)) * epsilon)); ++bl) {
  int block_index = to_sort[bl].second;
#else
  for (int block_index = 0; block_index < n_blocks; ++block_index) {
#endif

  int block_size = sosp.get_eigensystems()[block_index].eigenvalues.size();

  for (int state_index = 0; state_index < block_size; ++state_index) {
   state_t const& psi0 = sosp.get_eigensystems()[block_index].eigenstates[state_index];
   // do the first exp
   double dtau = (_begin == _end ? config->beta() : double(_begin->first));
   state_t psi = exp_h(psi0, dtau);
   
   for (auto it = _begin; it != _end;) { // do nothing if no operator
    // apply operator
    auto const& op = sosp.get_fundamental_operator_from_linear_index(it->second.dagger, it->second.linear_index);
    //auto const& op = sosp.get_fundamental_operator(it->second.dagger, it->second.block_index, it->second.inner_index);
    psi = op(psi);

    // apply exponential.
    double tau1 = double(it->first);
    ++it;
    dtau = (it == _end ? config->beta() : double(it->first)) - tau1;
    assert(dtau > 0);
    exp_h.apply(psi, dtau);
    // psi = exp_h (psi, dtau);
   }

   auto partial_trace = dot_product(psi0, psi);
   full_trace += partial_trace;
  }
 }
 return full_trace;
 }
}
