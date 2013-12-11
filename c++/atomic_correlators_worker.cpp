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
 auto _begin = config->oplist.rbegin();
 auto _end = config->oplist.rend();
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
  double tau1 = double(it->first);
  ++it;
  dtau = (it == _end ? config->beta() : double(it->first)) - tau1;
  assert(dtau > 0);

  bool one_non_zero = false;
  for (int n = 0; n < n_blocks; ++n) {
   if (blo[n] == -1) continue; // that chain has ended
   one_non_zero = true;
   blo[n] = sosp.fundamental_operator_connect(it1->second.dagger, it1->second.block_index, it1->second.inner_index, blo[n]);
   E_min_delta_tau[n] += dtau * sosp.get_eigensystems()[n].eigenvalues[0]; // delta_tau * E_min_of_the_block
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
#ifndef NO_FIRST_PASS
 for (int bl = 0; ((bl < n_bl) && (std::exp(- to_sort[bl].first) >= ( std::abs(full_trace)) * epsilon)); ++bl) {
  int block_index = to_sort[bl].second;
  //std::cout << "bl = " << bl << " exp (-Emin dtau) "<< std::exp(- to_sort[bl].first)  <<  " trace = "<< std::abs(full_trace)  << std::endl;
  //for (int block_index = n_blocks-1; block_index < n_blocks; ++block_index) {
#else
  for (int block_index = 0; block_index < n_blocks; ++block_index) {
#endif

  int block_size = sosp.get_eigensystems()[block_index].eigenvalues.size();

  for (int state_index = 0; state_index < block_size; ++state_index) {

   state_t const& psi0 = sosp.get_eigensystems()[block_index].eigenstates[state_index];

   // do the first exp
   double dtau = (_begin == _end ? config->beta() : double(_begin->first));
   state_t psi = exp_h(psi0, dtau);

   int cc = 0;
   bool psi_is_zero = false;
   for (auto it = _begin; it != _end; cc++) { // do nothing if no operator

    // apply operator
    auto const& op = sosp.get_fundamental_operator(it->second.dagger, it->second.block_index, it->second.inner_index);
    psi = op(psi);

    // std::cout << "gs energy interm "<< sosp.get_eigensystems()[psi.get_hilbert().get_index()].eigenvalues[0] << std::endl;
    // psi is already zero, makes no sense to proceed
    if (is_zero_state(psi)) {
     psi_is_zero = true;
     //if (cc !=0) std::cout << "Cancel after : "<< cc << std::endl ;
     break;
    }
    // if(is_zero_state(psi)) { if (cc !=0) std::cout << "Cancel after : "<< cc << std::endl ; return 0;}

    // apply exponential.
    double tau1 = double(it->first);
    ++it;
    dtau = (it == _end ? config->beta() : double(it->first)) - tau1;
    assert(dtau > 0);
    exp_h.apply(psi, dtau);
    // psi = exp_h (psi, dtau);
   }

   if (psi_is_zero) continue;

   auto partial_trace = dot_product(psi0, psi);
   full_trace += partial_trace;

   // WHY IS full_trace complex ??
   // std::cout  << "trace = "<< partial_trace << "   "<< full_trace << std::endl;

   // std::cout  << " partial trace "<< n << " "<< partial_trace<< std::endl;
   /*std::cout  << " number of BS "<< partial_traces.size() << std::endl;
     int i=0;
     for (auto x : partial_traces) if (std::abs(x)> 1.e-10) std::cout << i++ << "  "<< std::abs(x) << std::endl;
     std::cout  << "-----------"<<std::endl;
     */
  }
 }
 // std::cout << "trace = " << full_trace << std::endl;
 return full_trace;
 }
}
