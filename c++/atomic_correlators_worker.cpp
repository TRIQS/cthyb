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
     small_matrix_size(small_matrix_size) {
 histos.insert({"FirsTerm_FullTrace", {0, 10, 100, "hist_FirsTerm_FullTrace.dat"}});
 histos.insert({"FullTrace_ExpSumMin", {0, 10, 100, "hist_FullTrace_ExpSumMin.dat"}});
 histos.insert({"FullTrace_ExpSumMin", {0, 10, 100, "hist_FullTrace_ExpSumMin.dat"}});
 histo_bs_block = statistics::histogram{sosp.n_subspaces(), "hist_BS1.dat"};
}

//------------------------------------------------------------------------------

struct _p1 { double dtau; bool dag; long n;};

atomic_correlators_worker::result_t atomic_correlators_worker::operator()() {
 const auto _begin = config->oplist.crbegin();
 const auto _end = config->oplist.crend();
 const auto last_tau = config->beta();
 const int n_blocks = sosp.n_subspaces();
 const int config_size = config->oplist.size();

//#define NO_FIRST_PASS
#ifndef NO_FIRST_PASS
 // make a first pass to compute the bound for each term.
 std::vector<double> E_min_delta_tau(n_blocks, 0);
 std::vector<int> blo(n_blocks);

 std::vector<_p1> config_table(config_size);
 const double dtau0 = (_begin == _end ? config->beta() : double(_begin->first));

 // fill the table
 {
  int ii = 0;
  for (auto it = _begin; it != _end;) { // do nothing if no operator
   auto it1 = it;
   ++it;
   double dtau = (it == _end ? config->beta() : double(it->first)) - double(it1->first);
   config_table[ii++] = {dtau, it1->second.dagger, it1->second.linear_index};
  }
 }

 // first block only. MUST HAVE THE BLOCK IN Emin Order 
 // TEST ONLY
 double E_min_delta_tau0=0;

 bool TEST_PRI_ONLY = true;
 if (TEST_PRI_ONLY) {
  auto bl0 = 351;
  E_min_delta_tau0 = dtau0 * sosp.get_eigensystems()[351].eigenvalues[0];
  for (int i = 0; i < config_size; ++i) {
   // apply operator
   bl0 = sosp.fundamental_operator_connect_from_linear_index(config_table[i].dag, config_table[i].n, bl0);
   if (bl0 == -1) break;
   E_min_delta_tau0 += config_table[i].dtau  * sosp.get_eigensystems()[bl0].eigenvalues[0]; // delta_tau * E_min_of_the_block
  }
 }

 bool one_non_zero = false;
 for (int n = 0; n < n_blocks; ++n) {
  int bl = n;
  double sum_emin_dtau = dtau0 * sosp.get_eigensystems()[n].eigenvalues[0];
  for (int i = 0; i < config_size; ++i) {
   bl = sosp.fundamental_operator_connect_from_linear_index(config_table[i].dag, config_table[i].n, bl);
   if (bl == -1) break;
   sum_emin_dtau += config_table[i].dtau * sosp.get_eigensystems()[bl].eigenvalues[0]; // delta_tau * E_min_of_the_block
   if (TEST_PRI_ONLY && (sum_emin_dtau > E_min_delta_tau0 + 35)) {                                        // exp (-35) = 1.e-15
    bl = -1;
    break;
   }
  }
  blo[n] = bl;
  E_min_delta_tau[n] = sum_emin_dtau;
  one_non_zero |= (bl != -1);
 }
 if (!one_non_zero) return 0; // quick exit, the trace is structurally 0

 // Now sort the blocks
 std::vector<std::pair<double, int>> to_sort(n_blocks);
 int n_bl = 0; // the number of blocks giving non zero
 for (int n = 0; n < n_blocks; ++n)
  if (blo[n] == n)  // Must return to the SAME block, or trace is 0
   to_sort[n_bl++] = std::make_pair(E_min_delta_tau[n], n);

 std::sort(to_sort.begin(), to_sort.begin() + n_bl); // sort those vector

 // NOT much faster because the first part of the code IS LONGER
 // TOO MANY BLOCS ! --> do first GS, then search for better...
 //return std::exp(-to_sort[0].first); // QUICK estimate 

#endif

 result_t full_trace = 0;
 double epsilon = 1.e-15;
 double first_term = 0;

// To implement : regroup all the vector of the block for dgemm computation !
#ifndef NO_FIRST_PASS
 for (int bl = 0; ((bl < n_bl) && (std::exp(-to_sort[bl].first) >= (std::abs(full_trace)) * epsilon)); ++bl) {
  int block_index = to_sort[bl].second;
  auto exp_no_emin = std::exp(-to_sort[bl].first);
#else
 for (int block_index = 0; block_index < n_blocks; ++block_index) {
#endif

  int block_size = sosp.get_eigensystems()[block_index].eigenvalues.size();

  for (int state_index = 0; state_index < block_size; ++state_index) {
   state_t const& psi0 = sosp.get_eigensystems()[block_index].eigenstates[state_index];

   // do the first exp
   //dtau0 = (_begin == _end ? config->beta() : double(_begin->first));
   state_t psi = psi0;
   exp_h.apply_no_emin(psi, dtau0);

   for (auto it = _begin; it != _end;) { // do nothing if no operator
    // apply operator
    auto const& op = sosp.get_fundamental_operator_from_linear_index(it->second.dagger, it->second.linear_index);
    // auto const& op = sosp.get_fundamental_operator(it->second.dagger, it->second.block_index, it->second.inner_index);
    psi = op(psi);

    // apply exponential.
    double tau1 = double(it->first);
    ++it;
    double dtau = (it == _end ? config->beta() : double(it->first)) - tau1;
    assert(dtau > 0);
    exp_h.apply_no_emin(psi, dtau);
   }

   auto partial_trace_no_emin = dot_product(psi0, psi);
   auto partial_trace = partial_trace_no_emin * exp_no_emin;

   // CHECK conjecture
   if (std::abs(partial_trace_no_emin) > 1.0000001) throw "halte la !";

   /*
   if (bl==0) {
    std::cout << "-------" << std::endl;
    std::cout  << " block_index" << block_index<<std::endl;
    std::cout << "partial trace " << std::abs(partial_trace) << std::endl;
    std::cout << "partial_trace without emin" << std::abs(partial_trace_check) << std::endl;
    std::cout << "exp - sum emin" << std::exp(-to_sort[bl].first) << std::endl;
    std::cout << "exp - sum emin" << std::exp(-E_min_delta_tau[block_index]) << std::endl;
    std::cout << "<1 ?" << std::abs(partial_trace) / std::exp(-to_sort[bl].first) << std::endl;
   // std::cout << "partial_trace_noexp " << std::abs(partial_trace_noexp) << std::endl;
   }
 */

   if (bl == 0) first_term = partial_trace;
   full_trace += partial_trace;
  }
 }

 //bool use_histograms = false;
 bool use_histograms = true; //false;
 if (use_histograms) {
  auto abs_trace = std::abs(full_trace);
  if (abs_trace > 0) histos["FirsTerm_FullTrace"] << std::abs(first_term) / abs_trace;
  histos["FullTrace_ExpSumMin"] << std::abs(full_trace) / std::exp(-to_sort[0].first);
  histo_bs_block << to_sort[0].second;
 }

 return full_trace;
}
}
