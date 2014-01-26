#include "atomic_correlators_worker.hpp"
#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>
#include <algorithm>
#include <limits>

namespace cthyb_krylov {

atomic_correlators_worker::atomic_correlators_worker(configuration& c, sorted_spaces const& sosp_, double gs_energy_convergence,
                                                     int small_matrix_size, bool make_histograms, bool use_quick_trace_estimator,
                                                     int trace_estimator_n_blocks_guess, bool use_truncation, bool use_old_trace)
   : config(&c),
     sosp(sosp_),
     exp_h(sosp.get_hamiltonian(), sosp, gs_energy_convergence, small_matrix_size),
     small_matrix_size(small_matrix_size),
     make_histograms(make_histograms),
     use_quick_trace_estimator(use_quick_trace_estimator),
     trace_estimator_n_blocks_guess(trace_estimator_n_blocks_guess),
     use_truncation(use_truncation),
     use_old_trace(use_old_trace),
     time_spent_in_block(sosp.n_subspaces()),
     partial_over_full_trace(sosp.n_subspaces()),
     block_died_anal(sosp.n_subspaces(), 11),
     block_died_num(sosp.n_subspaces(), 11) {
 block_died_anal() = 0;
 block_died_num() = 0;
 if (make_histograms) {
  histos.insert({"FirsTerm_FullTrace", {0, 10, 100, "hist_FirsTerm_FullTrace.dat"}});
  histos.insert({"FullTrace_ExpSumMin", {0, 10, 100, "hist_FullTrace_ExpSumMin.dat"}});
  histos.insert({"FullTrace_over_Estimator", {0, 10, 100, "hist_FullTrace_over_Estimator.dat"}});
  histos.insert({"ExpBlock_over_ExpFirsTerm", {0, 1, 100, "hist_ExpBlock_over_ExpFirsTerm.dat"}});
  histo_bs_block = statistics::histogram{sosp.n_subspaces(), "hist_BS1.dat"};
  histo_trace_null_struc = statistics::histogram{2, "hist_trace_struct_nulle.dat"};

  // histo_opcount = statistics::histogram{100, "hist_opcount.dat"};

  for (int i = 0; i < 20; ++i) {
   std::stringstream s;
   s << "histo_n_blocks_after_steps" << i << ".dat";
   histo_n_blocks_after_steps.emplace_back(sosp.n_subspaces(), s.str());
  }

  for (int i = 0; i < 10; ++i) {
   std::stringstream s;
   s << "histo_opcount" << i << ".dat";
   histo_opcount.emplace_back(100, s.str());
  }
 }

 // insert the boundary points for the cache.
}

//------------------------------------------------------------------------------

atomic_correlators_worker::~atomic_correlators_worker() {

 boost::mpi::communicator world;
 std::string s = "time_and_partial_trace.dat";
 std::string t = "block_died_anal.dat"; // first x<=10 columns measure # times
 std::string u = "block_died_num.dat";  // block dies after x operators,
                                        // 11th column measures # times block
                                        // dies some point after 10 operators

 if (world.rank() == 0) {
  std::ofstream f(s);
  std::ofstream g(t);
  std::ofstream h(u);
  f << "Block  Time in block  Partial/Full Trace  Time*Partial/Full Trace" << std::endl;
  for (int i = 0; i < time_spent_in_block.size(); i++) {
   f << i << " " << time_spent_in_block[i] << " " << partial_over_full_trace[i] << " "
     << time_spent_in_block[i] * partial_over_full_trace[i] << std::endl;
   g << i;
   h << i;
   for (int j = 0; j < 11; j++) {
    g << " " << block_died_anal(i, j);
    h << " " << block_died_num(i, j);
   }
   g << std::endl;
   h << std::endl;
  }
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
 std::vector<_p1> config_table(config->size());
 const auto _begin = config->lowest_time_operator(); //config->oplist.crbegin();
 const auto _end = config->boundary_beta(); //config->oplist.crend();
 int ii = 0;
 for (auto it = _begin; it != _end;) { // do nothing if no operator
  auto it1 = it;
  --it;
  double dtau = double(it->first) - double(it1->first);
  config_table[ii++] = {dtau, it1->second.dagger, it1->second.linear_index};
 }
 return config_table;
}


// ------------------ trace estimate --------------

atomic_correlators_worker::result_t atomic_correlators_worker::estimate() {

 int config_size = config->size();
 auto dtau0 = double(config->lowest_time_operator()->first);
 double E_min_delta_tau_min = std::numeric_limits<double>::max();
 bool one_non_zero = false;
 auto config_table = make_config_table(config);
 int n_blocks = sosp.n_subspaces();
 int n_block_non_failed = 0;
 auto n_block_max = (trace_estimator_n_blocks_guess == -1 ? n_blocks : trace_estimator_n_blocks_guess);

 std::vector<int> n_blocks_after_steps(20, 0);

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
   if (i < 20) n_blocks_after_steps[i]++;
  }
  if (bl == n) {
   E_min_delta_tau_min = std::min(E_min_delta_tau_min, sum_emin_dtau);
   ++n_block_non_failed;
   one_non_zero = true;
  }
 } // loop over n

 // analysis
 // for (int i = 0; i < 10; ++i) std::cout << n_blocks_after_steps[i] << " "; std::cout << std::endl;
 // for (int i = 0; i < 10; ++i) histo_n_blocks_after_steps[i] << std::accumulate(n_blocks_after_steps.begin(),
 // n_blocks_after_steps.begin() + i + 1, 0);
 if (make_histograms)
  for (int i = 0; i < 20; ++i) histo_n_blocks_after_steps[i] << n_blocks_after_steps[i];

 if (!one_non_zero) return 0; // quick exit, the trace is structurally 0
 return std::exp(-E_min_delta_tau_min);
}

// ------------------- full trace computation -------------

atomic_correlators_worker::result_t atomic_correlators_worker::full_trace() {

 int n_blocks = sosp.n_subspaces();
 int config_size = config->size();
 std::vector<result_t> partial_trace_of_block(n_blocks);

 // make a first pass to compute the bound for each term.
 std::vector<double> E_min_delta_tau(n_blocks, 0);
 std::vector<int> blo(n_blocks);

 auto config_table = make_config_table(config);
 auto dtau0 = double(config->lowest_time_operator()->first);
 //double dtau0 = (_begin == _end ? config->beta() : double(_begin->first));

 std::vector<int> n_blocks_after_steps(20, 0);

 bool one_non_zero = false;
 double E_min_delta_tau_min = std::numeric_limits<double>::max() - 100;
 for (int n = 0; n < n_blocks; ++n) {
  int bl = n;
  double sum_emin_dtau = dtau0 * sosp.get_eigensystems()[n].eigenvalues[0];
  for (int i = 0; i < config_size; ++i) {
   bl = sosp.fundamental_operator_connect_from_linear_index(config_table[i].dag, config_table[i].n, bl);
   if (bl == -1) {
    block_died_anal(bl, std::max(i, 10))++;
    break;
   }
   sum_emin_dtau += config_table[i].dtau * sosp.get_eigensystems()[bl].eigenvalues[0]; // delta_tau * E_min_of_the_block
   if (sum_emin_dtau > E_min_delta_tau_min + 35) {                                     // exp (-35) = 1.e-15
    bl = -1;
    block_died_num(bl, std::max(i, 10))++;
    break;
   }
   if (i < 20) n_blocks_after_steps[i]++;
  }
  blo[n] = bl;
  E_min_delta_tau[n] = sum_emin_dtau;
  if (bl == n) {
   E_min_delta_tau_min = std::min(E_min_delta_tau_min, sum_emin_dtau);
   one_non_zero = true;
  }
 }

 if (make_histograms) histo_trace_null_struc << one_non_zero;

 if (!one_non_zero) return 0; // quick exit, the trace is structurally 0

 // analysis
 if (make_histograms) {
  std::vector<int> opcount(10, 0);
  for (auto const& p : *config) opcount[p.second.linear_index]++;
  for (int i = 0; i < 10; ++i) histo_opcount[i] << opcount[i] / 2;
  // histo_opcount << config_size/2; // histogram of the configuration size
  for (int i = 0; i < 20; ++i) histo_n_blocks_after_steps[i] << n_blocks_after_steps[i];

  if (0) {
   std::cout << config_size << " steps: ";
   for (int i = 0; i < 20; ++i) std::cout << n_blocks_after_steps[i] << " ";
   std::cout << std::endl; // <<*config << std::endl;
  }
 }
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

 auto exp_first_term = std::exp(-to_sort[0].first); // precompute largest term

 for (int bl = 0; ((bl < n_bl) && ((!use_truncation) || (std::exp(-to_sort[bl].first) >= (std::abs(full_trace)) * epsilon)));
      ++bl) {
  int block_index = to_sort[bl].second;
  auto exp_no_emin = std::exp(-to_sort[bl].first);

  // without the first pass would be
  // for (int block_index = 0; block_index < n_blocks; ++block_index)

  int block_size = sosp.get_eigensystems()[block_index].eigenvalues.size();

  if (make_histograms) {
   histos["ExpBlock_over_ExpFirsTerm"] << exp_no_emin / exp_first_term; // ratio of e^-E_bl to e^-E0
  }

  // -.-.-.-.-.-.-.-.-.   Old implementation of the trace -.-.-.-.-.-.-
  if (use_old_trace) {

  for (int state_index = 0; state_index < block_size; ++state_index) {
    state_t const& psi0 = sosp.get_eigensystems()[block_index].eigenstates[state_index];

    // do the first exp
    // dtau0 = (_begin == _end ? config->beta() : double(_begin->first));
    state_t psi = psi0;
    exp_h.apply_no_emin(psi, dtau0);

    auto _begin = config->lowest_time_operator();
    auto _end = config->boundary_zero();
    for (auto it = _begin; it != _end;) { // do nothing if no operator
     // apply operator
     auto const& op = sosp.get_fundamental_operator_from_linear_index(it->second.dagger, it->second.linear_index);
     psi = op(psi);

     // apply exponential.
     double tau1 = double(it->first);
     --it;
     double dtau = double(it->first) - tau1;
     assert(dtau > 0);
     exp_h.apply_no_emin(psi, dtau);
    }

    auto partial_trace_no_emin = dot_product(psi0, psi);
    auto partial_trace = partial_trace_no_emin * exp_no_emin;
    if (std::abs(partial_trace_no_emin) > 1.0000001) throw "halte la !"; // CHECK conjecture

    if (bl == 0) first_term = partial_trace;
    full_trace += partial_trace;
   }
  } else { // -.-.-.-.-.-.-.-.-.  New implementation of the trace -.-.-.-.-.-.-

   auto const& eigensystem = sosp.get_eigensystems()[block_index];
   auto space_dim = eigensystem.eigenvalues.size();
   auto const& U_boundary = sosp.get_eigensystems()[block_index].unitary_matrix;
   auto M = U_boundary;
   // auto M2 = triqs::arrays::matrix<double>(100,100);//  CHANGE THIS !!
   auto _ = range{};
   // do the first exponential
   for (int n = 1; n < space_dim; ++n) M(n, _) *= exp(-dtau0 * (eigensystem.eigenvalues(n) - eigensystem.eigenvalues(0)));

   auto B = block_index;

   for (int i = 0; i < config_size; ++i) {
    // get the next block
    auto Bp = sosp.fundamental_operator_connect_from_linear_index(config_table[i].dag, config_table[i].n, B);
    if (Bp == -1) TRIQS_RUNTIME_ERROR << " Internal error : Bp =-1";
    // apply operator
    auto const& Cop = sosp.fundamental_operator_matrix_from_linear_index(config_table[i].dag, config_table[i].n, B);
    auto const& eigensystem = sosp.get_eigensystems()[Bp];
    auto space_dim_p = eigensystem.eigenvalues.size();

    if ((space_dim == 1) && (space_dim_p == 1)) // some quick optimisation
     M *= Cop(0, 0);
    else {
     auto M2 = Cop * M;
     swap(M2, M);
     // improve this with a larger matrix to avoid allocation ???
     // triqs::arrays::blas::gemm(1.0, Cop, M, 0.0, M2(range(0,space_dim_p), range(0,space_dim));
    }
    B = Bp;
    space_dim = space_dim_p;

    // apply exponential.
    if (space_dim_p == 1) continue;
    for (int n = 1; n < space_dim_p; ++n)
     M(n, _) *= std::exp(-config_table[i].dtau * (eigensystem.eigenvalues(n) - eigensystem.eigenvalues(0)));

    // measure time spent in block B
    time_spent_in_block[B] += double(config_table[i].dtau);

   } // loop on the c ops of the configuration

   auto partial_trace_no_emin = trace(U_boundary.transpose() * M);
   auto partial_trace = partial_trace_no_emin * exp_no_emin;
   if (std::abs(partial_trace_no_emin) > space_dim * 1.0000001) throw "halte la !"; // CHECK conjecture

   if (bl == 0) first_term = partial_trace;
   full_trace += partial_trace;
   partial_trace_of_block[block_index] = partial_trace;

  } // -.-.-.-.-.-.-.-.-.  choice of trace computation method -.-.-.-.-.-.-

 } // end of loop on blocks

 if (make_histograms) {
  auto abs_trace = std::abs(full_trace);
  if (abs_trace > 0) histos["FirsTerm_FullTrace"] << std::abs(first_term) / abs_trace;
  histos["FullTrace_ExpSumMin"] << std::abs(full_trace) / std::exp(-to_sort[0].first);
  histo_bs_block << to_sort[0].second;
 }
 for (int i = 0; i < partial_trace_of_block.size(); i++) partial_over_full_trace[i] += partial_trace_of_block[i] / full_trace;
 return full_trace;
}
}
