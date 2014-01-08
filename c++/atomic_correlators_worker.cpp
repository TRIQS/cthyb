#include "atomic_correlators_worker.hpp"
#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>
#include <algorithm>
#include <limits>

namespace cthyb_krylov {

atomic_correlators_worker::atomic_correlators_worker(configuration& c, sorted_spaces const& sosp_, double gs_energy_convergence,
                                                     int small_matrix_size, bool make_histograms, bool use_quick_trace_estimator,
                                                     int trace_estimator_n_blocks_guess)
   : config(&c),
     sosp(sosp_),
     exp_h(sosp.get_hamiltonian(), sosp, gs_energy_convergence, small_matrix_size),
     small_matrix_size(small_matrix_size),
     make_histograms(make_histograms),
     use_quick_trace_estimator(use_quick_trace_estimator),
     trace_estimator_n_blocks_guess(trace_estimator_n_blocks_guess) {
 if (make_histograms) {
  histos.insert({"FirsTerm_FullTrace", {0, 10, 100, "hist_FirsTerm_FullTrace.dat"}});
  histos.insert({"FullTrace_ExpSumMin", {0, 10, 100, "hist_FullTrace_ExpSumMin.dat"}});
  histos.insert({"FullTrace_over_Estimator", {0, 10, 100, "hist_FullTrace_over_Estimator.dat"}});
  histo_bs_block = statistics::histogram{sosp.n_subspaces(), "hist_BS1.dat"};
  histo_block_size = statistics::histogram{100, "hist_block_size.dat"};
 }
}

//------------------------------------------------------------------------------


atomic_correlators_worker::result_t atomic_correlators_worker::operator()() {
 return (use_quick_trace_estimator ? estimate() : full_trace());
}

//------------------------------------------------------------------------------

atomic_correlators_worker::result_t atomic_correlators_worker::full_trace_over_estimator() {
 if (!use_quick_trace_estimator) return 1;
 auto r = full_trace() / estimate();
 if (make_histograms) histos["FullTrace_over_Estimator"] << std::abs(r);
 return r;
}

//------------------------- make_config_table ------------------------------------

struct _p1 {
 double dtau;
 bool dag;
 long n;
};

// a small temporary table from the config, to optimize the sweep in the algorithms below
std::vector<_p1> make_config_table(const configuration* config) {
 std::vector<_p1> config_table(config->oplist.size());
 const auto _begin = config->oplist.crbegin();
 const auto _end = config->oplist.crend();
 int ii = 0;
 for (auto it = _begin; it != _end;) { // do nothing if no operator
  auto it1 = it;
  ++it;
  double dtau = (it == _end ? config->beta() : double(it->first)) - double(it1->first);
  config_table[ii++] = {dtau, it1->second.dagger, it1->second.linear_index};
 }
 return config_table;
}

// ------------------ trace estimate --------------

atomic_correlators_worker::result_t atomic_correlators_worker::estimate() {

 auto _begin = config->oplist.crbegin();
 auto _end = config->oplist.crend();
 int config_size = config->oplist.size();
 double dtau0 = (_begin == _end ? config->beta() : double(_begin->first));
 double E_min_delta_tau_min = std::numeric_limits<double>::max();
 bool one_non_zero = false;
 auto config_table = make_config_table(config);
 int n_blocks = sosp.n_subspaces();
 int n_block_non_failed = 0;
 auto n_block_max = (trace_estimator_n_blocks_guess == -1 ? n_blocks : trace_estimator_n_blocks_guess);

 for (int n = 0; (n < n_blocks) && (n_block_non_failed < n_block_max); ++n) {
  int bl = n;
  double sum_emin_dtau = dtau0 * sosp.get_eigensystems()[n].eigenvalues[0];
  for (int i = 0; i < config_size; ++i) {
   bl = sosp.fundamental_operator_connect_from_linear_index(config_table[i].dag, config_table[i].n, bl);
   if (bl == -1) break;
   sum_emin_dtau += config_table[i].dtau * sosp.get_eigensystems()[bl].eigenvalues[0]; // delta_tau * E_min_of_the_block
   if (sum_emin_dtau > E_min_delta_tau_min) {
    bl = -1;
    break;
   }
  }
  if (bl != -1) {
   E_min_delta_tau_min = std::min(E_min_delta_tau_min, sum_emin_dtau);
   ++n_block_non_failed;
   one_non_zero = true;
  }
 }                            // loop over n
 if (!one_non_zero) return 0; // quick exit, the trace is structurally 0
 return std::exp(-E_min_delta_tau_min);
}

// ------------------- full trace computation -------------

atomic_correlators_worker::result_t atomic_correlators_worker::full_trace() {

 auto _begin = config->oplist.crbegin();
 auto _end = config->oplist.crend();
 int n_blocks = sosp.n_subspaces();
 int config_size = config->oplist.size();

 // make a first pass to compute the bound for each term.
 std::vector<double> E_min_delta_tau(n_blocks, 0);
 std::vector<int> blo(n_blocks);

 auto config_table = make_config_table(config);
 double dtau0 = (_begin == _end ? config->beta() : double(_begin->first));

 bool one_non_zero = false;
 double E_min_delta_tau_min = std::numeric_limits<double>::max() - 100;
 for (int n = 0; n < n_blocks; ++n) {
  int bl = n;
  double sum_emin_dtau = dtau0 * sosp.get_eigensystems()[n].eigenvalues[0];
  for (int i = 0; i < config_size; ++i) {
   bl = sosp.fundamental_operator_connect_from_linear_index(config_table[i].dag, config_table[i].n, bl);
   if (bl == -1) break;
   sum_emin_dtau += config_table[i].dtau * sosp.get_eigensystems()[bl].eigenvalues[0]; // delta_tau * E_min_of_the_block
   if (sum_emin_dtau > E_min_delta_tau_min + 35) {                                     // exp (-35) = 1.e-15
    bl = -1;
    break;
   }
  }
  blo[n] = bl;
  E_min_delta_tau[n] = sum_emin_dtau;
  if (bl != -1) E_min_delta_tau_min = std::min(E_min_delta_tau_min, sum_emin_dtau);
  one_non_zero |= (bl != -1);
 }
 if (!one_non_zero) return 0; // quick exit, the trace is structurally 0

 // Now sort the blocks
 std::vector<std::pair<double, int>> to_sort(n_blocks);
 int n_bl = 0; // the number of blocks giving non zero
 for (int n = 0; n < n_blocks; ++n)
  if (blo[n] == n) // Must return to the SAME block, or trace is 0
   to_sort[n_bl++] = std::make_pair(E_min_delta_tau[n], n);

 std::sort(to_sort.begin(), to_sort.begin() + n_bl); // sort those vector

 // - end first pass

 result_t full_trace = 0;
 double epsilon = 1.e-15;
 double first_term = 0;

 // To implement : regroup all the vector of the block for dgemm computation !
 for (int bl = 0; ((bl < n_bl) && (std::exp(-to_sort[bl].first) >= (std::abs(full_trace)) * epsilon)); ++bl) {
  int block_index = to_sort[bl].second;
  auto exp_no_emin = std::exp(-to_sort[bl].first);

  // without the first pass would be
  // for (int block_index = 0; block_index < n_blocks; ++block_index)

  int block_size = sosp.get_eigensystems()[block_index].eigenvalues.size();

  if (make_histograms) histo_block_size << block_size;

  for (int state_index = 0; state_index < block_size; ++state_index) {
   state_t const& psi0 = sosp.get_eigensystems()[block_index].eigenstates[state_index];

   // do the first exp
   // dtau0 = (_begin == _end ? config->beta() : double(_begin->first));
   state_t psi = psi0;
   exp_h.apply_no_emin(psi, dtau0);

   for (auto it = _begin; it != _end;) { // do nothing if no operator
    // apply operator
    auto const& op = sosp.get_fundamental_operator_from_linear_index(it->second.dagger, it->second.linear_index);
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
   if (std::abs(partial_trace_no_emin) > 1.0000001) throw "halte la !"; // CHECK conjecture

   if (bl == 0) first_term = partial_trace;
   full_trace += partial_trace;
  }
 }
 if (make_histograms) {
  auto abs_trace = std::abs(full_trace);
  if (abs_trace > 0) histos["FirsTerm_FullTrace"] << std::abs(first_term) / abs_trace;
  histos["FullTrace_ExpSumMin"] << std::abs(full_trace) / std::exp(-to_sort[0].first);
  histo_bs_block << to_sort[0].second;
 }
 return full_trace;
}
}
