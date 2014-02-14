#include "atomic_correlators_worker.hpp"
#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>
#include <algorithm>
#include <limits>

namespace cthyb_krylov {

atomic_correlators_worker::atomic_correlators_worker(configuration& c, sorted_spaces const& sosp_, utility::parameters const& p)
   : config(&c),
     sosp(sosp_),
     exp_h(sosp.get_hamiltonian(), sosp, p["krylov_gs_energy_convergence"], p["krylov_small_matrix_size"]),
     time_spent_in_block(sosp.n_subspaces(),0),
     partial_over_full_trace(sosp.n_subspaces(),0) {

 make_histograms = p["make_path_histograms"];
 use_truncation = p["use_truncation"];
 use_old_trace = p["use_old_trace"];

 std::string ms = p["trace_estimator"];
 //std::string ms = utility::extract<std::string>(p["trace_estimator"]);
 try {
  estimator_method = std::map<std::string, estimator_method_t>{{"None", estimator_method_t::None},
                                                               {"Simple", estimator_method_t::Simple},
                                                               {"WithCache", estimator_method_t::WithCache},
                                                               {"TraceEpsilon", estimator_method_t::TraceEpsilon},
                                                               {"Experimental1", estimator_method_t::Experimental1}}.at(ms);
 }
 catch (...) {
  TRIQS_RUNTIME_ERROR << "trace_estimator method " << ms << "not recognized";
 }
 
 //int max_subspace_dim=0;
 //for (int nsp = 0; nsp < sosp.n_subspaces(); ++nsp) max_subspace_dim = std::max(max_subspace_dim, sosp.subspace(nsp).dimension());
 //M_work.resize(max_subspace_dim, max_subspace_dim);
 //M_work2.resize(max_subspace_dim, max_subspace_dim);

 if (make_histograms) {
  histos.insert({"FirsTerm_FullTrace", {0, 10, 100, "hist_FirsTerm_FullTrace.dat"}});
  histos.insert({"FullTrace_ExpSumMin", {0, 10, 100, "hist_FullTrace_ExpSumMin.dat"}});
  histos.insert({"FullTrace_over_Estimator", {0, 10, 100, "hist_FullTrace_over_Estimator.dat"}});
  histos.insert({"ExpBlock_over_ExpFirsTerm", {0, 1, 100, "hist_ExpBlock_over_ExpFirsTerm.dat"}});
  histo_bs_block = statistics::histogram{sosp.n_subspaces(), "hist_BS1.dat"};
  histo_trace_null_struc = statistics::histogram{2, "hist_trace_struct_nulle.dat"};
  histo_n_block_kept = statistics::histogram{1000, "histo_n_block_kept.dat"};
  histo_n_block_at_end = statistics::histogram{1000, "histo_n_block_at_end.dat"};
  histo_n_block_after_esti = statistics::histogram{1000, "histo_n_block_after_esti.dat"};
  histo_blocks_after_esti = statistics::histogram{sosp.n_subspaces(), "histo_blocks_after_esti.dat"};
 
  for (int i = 0; i < 50; ++i) {
   std::stringstream s;
   s << "histo_n_blocks_after_steps" << i << ".dat";
   histo_n_blocks_after_steps.emplace_back(sosp.n_subspaces(), s.str());
  }

  for (int i = 0; i < 50; ++i) {
   std::stringstream s;
   s << "histo_n_blocks_cache_rl_" << i << ".dat";
   histo_n_blocks_cache_rl.emplace_back(sosp.n_subspaces(), s.str());
  }
  for (int i = 0; i < 50; ++i) {
   std::stringstream s;
   s << "histo_n_blocks_cache_lr_" << i << ".dat";
   histo_n_blocks_cache_lr.emplace_back(sosp.n_subspaces(), s.str());
  }

  for (int i = 0; i < n_orbitals; ++i) {
   std::stringstream s;
   s << "histo_opcount" << i << ".dat";
   histo_opcount.emplace_back(100, s.str());
  }

  trunc_block = std::vector<int>{10,20,30,40,50,100};
  for (int i : trunc_block) {
   std::stringstream s;
   std::stringstream t;
   s << "TruncatedTrace_over_FullTrace" << i;
   t << "hist_TruncatedTrace_over_FullTrace" << i << ".dat";
   histos.insert({s.str(), {0, 3, 300, t.str()}});
  }
 }

 // must be AFTER the init of the histogram ! 
 cache_update();
}

//------------------------------------------------------------------------------

atomic_correlators_worker::~atomic_correlators_worker() {

 if (make_histograms) {
  boost::mpi::communicator world;
  std::string s = "hist_time_and_partial_trace.dat";
 
  if (world.rank() == 0) {
   std::ofstream f(s);
   f << "Block  Time in block  Partial/Full Trace  Time*Partial/Full Trace" << std::endl;
   for (int i = 0; i < time_spent_in_block.size(); i++) {
    f << i << " " << time_spent_in_block[i] << " " << partial_over_full_trace[i] << " "
      << time_spent_in_block[i] * partial_over_full_trace[i] << std::endl;
   }
  }
 }
}

//------------------------------------------------------------------------------

atomic_correlators_worker::result_t atomic_correlators_worker::estimate(time_pt t1, time_pt t2) {

 // DEBUG ONLY
 //auto r1 = full_trace(); 
 //auto r2 = full_trace2(); 
 ////if (std::abs(r1/r2 -1) > 1.e-10) 
 //if (std::abs(r1-r2)> (1.e-10 + 0.1 * std::abs(r1)) ) {
 //std::cout << r1 << " == " << r2 << " |r1 -r2| "<< std::abs(r1-r2) << std::endl;
 // throw "ERROR";
 // }
 //return r1;

 switch (estimator_method) {
  case estimator_method_t::None:
   return full_trace();
  case estimator_method_t::TraceEpsilon:
   return full_trace(trace_epsilon_estimator());
  case estimator_method_t::Experimental1:
   return full_trace2();
  case estimator_method_t::Simple:
   return estimate_simple();
  case estimator_method_t::WithCache:
   return estimate_with_cache(t1, t2);
 }
}

//------------------------------------------------------------------------------


atomic_correlators_worker::result_t atomic_correlators_worker::estimate() {

 if (estimator_method == estimator_method_t::None) return full_trace();
 if (estimator_method == estimator_method_t::TraceEpsilon) return full_trace(trace_epsilon_estimator());
 if (estimator_method == estimator_method_t::Experimental1) return full_trace2();
 return estimate_simple();
}

//------------------------------------------------------------------------------

atomic_correlators_worker::result_t atomic_correlators_worker::full_trace_over_estimator() {
 
 atomic_correlators_worker::result_t r=1;
 
 switch (estimator_method) {
  case estimator_method_t::None:
  case estimator_method_t::Experimental1:
   break;
  case estimator_method_t::TraceEpsilon:
   r = full_trace() / full_trace(trace_epsilon_estimator());
   break;
  case estimator_method_t::Simple:
  case estimator_method_t::WithCache:
   r = full_trace() / estimate_simple();
 }

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
 const auto _begin = config->lowest_time_operator(); // config->oplist.crbegin();
 const auto _end = config->boundary_beta();          // config->oplist.crend();
 int ii = 0;
 for (auto it = _begin; it != _end;) { // do nothing if no operator
  auto it1 = it;
  --it;
  double dtau = double(it->first) - double(it1->first);
  config_table[ii++] = {dtau, it1->second.dagger, it1->second.linear_index};
 }
 return config_table;
}

// ------------------ trace estimate with cache --------------

// tau1, tau2 : time of the last changes since last cache update.
// the trace is recomputed by windowing tau1, tau2 with 2 operators
// then gluing together the r, l cached computations of these operators
// and then running the estimate computation in the middle,
// breaking in particular when above the current maximum.
// todo : accelerate the iteration over the config with a table ?
//
atomic_correlators_worker::result_t atomic_correlators_worker::estimate_with_cache(time_pt tau1, time_pt tau2) {

 auto tl = std::max(tau1, tau2);
 auto tr = std::min(tau1, tau2);
 
 // The operators just at the left (higher time) than tl
 auto opl = config->operator_just_after(tl);
 auto opr = config->operator_just_before(tr);
 
 auto& c_l = cache.at(opl->first).l;
 auto& c_r = cache.at(opr->first).r;

 double E_min_delta_tau_min = std::numeric_limits<double>::max();
 bool one_non_zero = false;
 const int n_blocks = sosp.n_subspaces();
 int n_block_kept = 0;
 int n_block_at_end = 0;

 for (int n = 0; n < n_blocks; ++n) {
  if ((c_l[n].current_block_number == -1) || (c_r[n].current_block_number == -1)) continue;
  n_block_kept++;
  double sum_emin_dtau = c_l[n].emin_dtau_acc + c_r[n].emin_dtau_acc;
  int bl = c_l[n].current_block_number;
  auto op = opl;
  ++op;
  double tp = double(op->first);
  sum_emin_dtau += (double(opl->first) - tp) * sosp.get_eigensystems()[bl].eigenvalues[0]; // delta_tau * E_min_of_the_block
  while (op != opr) {
   bl = sosp.fundamental_operator_connect_from_linear_index(!op->second.dagger, op->second.linear_index, bl);
   if (bl == -1) break;
   ++op;
   sum_emin_dtau += (tp - op->first) * sosp.get_eigensystems()[bl].eigenvalues[0]; // delta_tau * E_min_of_the_block
   if (sum_emin_dtau > E_min_delta_tau_min) {
    bl = -1;
    break;
   }
   tp = double(op->first);
  }
  if ((bl != -1) && (bl == c_r[n].current_block_number)) {
   n_block_at_end++;
   E_min_delta_tau_min = std::min(E_min_delta_tau_min, sum_emin_dtau);
   one_non_zero = true;
  }
 } // loop over n

 histo_n_block_kept << n_block_kept;
 histo_n_block_at_end << n_block_at_end;
 if (!one_non_zero) return 0; // the trace is structurally 0
 //if (std::abs(E_min_delta_tau_min - estimate_simple(true) )>1.e-5) TRIQS_RUNTIME_ERROR << " FATAL";
 return std::exp(-E_min_delta_tau_min);
}

// ------------------ refresh the cache completely : only done when accepted, which is rare --------------

void atomic_correlators_worker::cache_update() {
 // only if the estimator method is with cache do we do something...
 if (estimator_method != estimator_method_t::WithCache) return; 
 cache_update_impl();
}
 
void atomic_correlators_worker::cache_update_impl() {

 // a brutal solution as first implementation : clean the cache and rebuild
 // better to add the cache in the config ?
 cache.clear();
 // insert the boundary points for the cache.
 {
  auto c = make_cache_point();
  for (int n = 0; n < sosp.n_subspaces(); ++n) {
   c.r[n].current_block_number = n;
   c.l[n].current_block_number = n;
  }
  cache.insert({time_pt::make_zero(config->beta()), c});
  cache.insert({time_pt::make_beta(config->beta()), c});
 }

 for (auto const& op : *config) cache.insert({op.first, cache_point_t{sosp.n_subspaces()}});

 std::vector<int> n_blocks_lr(histo_n_blocks_cache_lr.size(), 0);
 std::vector<int> n_blocks_rl(histo_n_blocks_cache_rl.size(), 0);
 
 const int n_blocks = sosp.n_subspaces();
 for (int n = 0; n < n_blocks; ++n) {

  // compute from beta -> 0
  int bl = n, i = 0;
  double emin_dtau_acc = 0;
  for (auto l = config->boundary_beta(), r = config->highest_time_operator(), _rend = config->boundary_zero(); r != _rend;
       ++r, ++l) {
   // evolve and act with the operator from beta -> 0
   auto& c = cache.at(r->first).l;
   emin_dtau_acc += (l->first - r->first) * sosp.get_eigensystems()[bl].eigenvalues[0]; // delta_tau * E_min_of_the_block
   bl = sosp.fundamental_operator_connect_from_linear_index(!r->second.dagger, r->second.linear_index, bl);
   c[n].current_block_number = bl;
   c[n].emin_dtau_acc = emin_dtau_acc;
   if (bl == -1) break;
   // instrumentation
   ++i;
   if (i < n_blocks_lr.size()) n_blocks_lr[i]++;
  }

  // compute from 0 -> beta
  bl = n;
  i = 0;
  emin_dtau_acc = 0;
  for (auto l = config->lowest_time_operator(), r = config->boundary_zero(), _lend = config->boundary_beta(); l != _lend;
       --r, --l) {
   // evolve and act with the operator from 0 -> beta
   auto& c = cache.at(l->first).r;
   emin_dtau_acc += (l->first - r->first) * sosp.get_eigensystems()[bl].eigenvalues[0]; // delta_tau * E_min_of_the_block
   bl = sosp.fundamental_operator_connect_from_linear_index(l->second.dagger, l->second.linear_index, bl);
   c[n].current_block_number = bl;
   c[n].emin_dtau_acc = emin_dtau_acc;
   if (bl == -1) break;
   // instrumentation
   ++i;
   if (i < n_blocks_rl.size()) n_blocks_rl[i]++;
   }
 }

 // analysis
 if (make_histograms) {
  std::vector<int> opcount(n_orbitals, 0);  // maximum number of orbitals is n_orbitals
  for (auto const& p : *config) opcount[p.second.linear_index]++;
  for (int i = 0; i < n_orbitals; ++i) histo_opcount[i] << opcount[i] / 2;
  // histo_opcount << config_size/2; // histogram of the configuration size
  for (int i = 0; i < std::min(config->size(),int(histo_n_blocks_cache_lr.size())); ++i) {
   //std::cout << i << histo_n_blocks_cache_lr.size() << " " << histo_n_blocks_cache_rl.size() << std::endl;
   histo_n_blocks_cache_lr[i] << n_blocks_lr[i];
   histo_n_blocks_cache_rl[i] << n_blocks_rl[i];
  }
 }

}

// ------------------ trace estimate --------------

atomic_correlators_worker::result_t atomic_correlators_worker::estimate_simple(bool no_exp) {

 int config_size = config->size();
 auto dtau0 = double(config->lowest_time_operator()->first);
 double E_min_delta_tau_min = std::numeric_limits<double>::max();
 bool one_non_zero = false;
 auto config_table = make_config_table(config);
 int n_blocks = sosp.n_subspaces();

 std::vector<int> n_blocks_after_steps(50, 0);

 for (int n = 0; (n < n_blocks); ++n) {
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
   if (i < 50) n_blocks_after_steps[i]++;
  }
  if (bl == n) {
   E_min_delta_tau_min = std::min(E_min_delta_tau_min, sum_emin_dtau);
   one_non_zero = true;
  }
 } // loop over n

 // analysis
 // for (int i = 0; i < 10; ++i) std::cout << n_blocks_after_steps[i] << " "; std::cout << std::endl;
 // for (int i = 0; i < 10; ++i) histo_n_blocks_after_steps[i] << std::accumulate(n_blocks_after_steps.begin(),
 // n_blocks_after_steps.begin() + i + 1, 0);
 if (make_histograms)
  for (int i = 0; i < 50; ++i) histo_n_blocks_after_steps[i] << n_blocks_after_steps[i];

 if (!one_non_zero) return 0; // quick exit, the trace is structurally 0
 
 if (no_exp) return E_min_delta_tau_min;
 return std::exp(-E_min_delta_tau_min);
}

// ------------

// TO BE MOVE TO THE TRIQS LIB : norm1 of a matrix
template <typename MatrixType> typename MatrixType::value_type norm1(MatrixType const& m) {
 auto res = typename MatrixType::value_type{};
 foreach(m, [&res, &m](int i, int j) {
  using std::abs;
  res += abs(m(i, j));
 });
 return res;
}

// ------------------- full trace computation -------------

atomic_correlators_worker::result_t atomic_correlators_worker::full_trace(double epsilon) {

 int n_blocks = sosp.n_subspaces();
 int config_size = config->size();
 double log_epsilon = -std::log(epsilon);
 auto config_table = make_config_table(config);
 auto dtau0 = double(config->lowest_time_operator()->first);
 
 // arrays for histo
 std::vector<result_t> partial_trace_of_block(n_blocks,0);
 std::vector<result_t> partial_trace_up_to_block(trunc_block.size(),0);
 std::vector<int> n_blocks_after_steps(50, 0); 
 
 // make a first pass to compute the bound for each term.
 std::vector<double> E_min_delta_tau(n_blocks, 0);
 std::vector<bool> is_block_kept(n_blocks, false);
 bool one_non_zero = false;
 double E_min_delta_tau_min = std::numeric_limits<double>::max() - 100;

for ( int uu=0; uu<1; ++uu) { // JSUT A TRICK TO EVALUATE TEH SPEED OF THIS : put uu < 2 or 3 
 one_non_zero = false;
 E_min_delta_tau_min = std::numeric_limits<double>::max() - 100;
 
 for (int n = 0; n < n_blocks; ++n) {
  int bl = n;
  double sum_emin_dtau = dtau0 * sosp.get_eigensystems()[n].eigenvalues[0];
  for (int i = 0; i < config_size; ++i) {
   bl = sosp.fundamental_operator_connect_from_linear_index(config_table[i].dag, config_table[i].n, bl);
   if (bl == -1) break;
   sum_emin_dtau += config_table[i].dtau * sosp.get_eigensystems()[bl].eigenvalues[0]; // delta_tau * E_min_of_the_block
   if (use_truncation && (sum_emin_dtau > E_min_delta_tau_min + log_epsilon)) {        // exp (-35) = 1.e-15
    bl = -1;
    break;
    }
   if (make_histograms && (i < 50)) n_blocks_after_steps[i]++;
  }
  E_min_delta_tau[n] = sum_emin_dtau;
  if (bl == n) {// Must return to the SAME block, or trace is 0
   is_block_kept[n] = true;
   E_min_delta_tau_min = std::min(E_min_delta_tau_min, sum_emin_dtau);
   one_non_zero = true;
  }
 }
}

 if (make_histograms) histo_trace_null_struc << one_non_zero;

 if (!one_non_zero) return 0; // quick exit, the trace is structurally 0

 // Now sort the blocks
 std::vector<std::pair<double, int>> to_sort(n_blocks);
 int n_bl = 0; // the number of blocks giving non zero
 for (int n = 0; n < n_blocks; ++n)
  if (is_block_kept[n] && (!use_truncation || (E_min_delta_tau[n] < E_min_delta_tau_min + log_epsilon))) // cut if too small
   to_sort[n_bl++] = std::make_pair(E_min_delta_tau[n], n);

 std::sort(to_sort.begin(), to_sort.begin() + n_bl);

 // analysis
 if (make_histograms) {
  std::vector<int> opcount(n_orbitals, 0);  // maximum number of orbitals is n_orbitals
  for (auto const& p : *config) opcount[p.second.linear_index]++;
  for (int i = 0; i < n_orbitals; ++i) histo_opcount[i] << opcount[i] / 2;
  // histo_opcount << config_size/2; // histogram of the configuration size

  for (int i = 0; i < std::min(config->size(), int(histo_n_blocks_after_steps.size())); ++i)
   histo_n_blocks_after_steps[i] << n_blocks_after_steps[i];

  // DEBUG
  //cache_update_impl();
  
  histo_n_block_after_esti << n_bl;// number of blocks after the estimator
  for (auto const & x: to_sort) histo_blocks_after_esti << x.second;// the blocks after the estimator

  //for (auto const& x : to_sort)
  // if (x.second != 0) std::cout << " block : " << x.second << " emin_dtau " << x.first << " min "<< to_sort[0].first<<  std::endl;
 }

 // ------------------- end first pass

 result_t full_trace = 0;
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
   auto space_dim_boundary = eigensystem.eigenvalues.size();
   auto space_dim = space_dim_boundary;
   auto const& U_boundary = sosp.get_eigensystems()[block_index].unitary_matrix;
   auto M = U_boundary;
   //M_work(range(0,space_dim), range(0,space_dim)) = U_boundary;
   auto _ = range{};
   // do the first exponential
   for (int n = 1; n < space_dim; ++n) M(n, _) *= exp(-dtau0 * (eigensystem.eigenvalues(n) - eigensystem.eigenvalues(0)));
   //for (int n = 1; n < space_dim; ++n) M_work(n, _) *= exp(-dtau0 * (eigensystem.eigenvalues(n) - eigensystem.eigenvalues(0)));

   auto B = block_index;

   bool loop_broken = false;
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
     //M_work(...) *= Cop(0, 0);
    else {
     //triqs::arrays::blas::gemm(1.0, Cop, M_work(range(0,space_dim), range(0,space_dim_boundary)), 0.0, M_work2(range(0,space_dim_p), range(0,space_dim_boundary)));
     auto M2 = Cop * M;
     swap(M, M2);
     //swap(M_work, M_work2);
     // improve this with a larger matrix to avoid allocation ???
     // triqs::arrays::blas::gemm(1.0, Cop, M, 0.0, M2(range(0,space_dim_p), range(0,space_dim));
    }
    B = Bp;
    space_dim = space_dim_p;

    // apply exponential.
    if (space_dim_p != 1) {
     for (int n = 1; n < space_dim_p; ++n)
      M(n, _) *= std::exp(-config_table[i].dtau * (eigensystem.eigenvalues(n) - eigensystem.eigenvalues(0)));
      //M_work(n, _) *= std::exp(-config_table[i].dtau * (eigensystem.eigenvalues(n) - eigensystem.eigenvalues(0)));
    }
    // measure time spent in block B
    if (make_histograms) time_spent_in_block[B] += double(config_table[i].dtau);

    // further cut : if what we are computing is in fact already smaller than epsilon * full_trace
    //if (use_truncation && (norm1(M_work(range(0,space_dim),range(0,space_dim))) * exp_no_emin < epsilon * std::abs(full_trace))) {
    if (use_truncation && (norm1(M) * exp_no_emin < epsilon * std::abs(full_trace))) {
     loop_broken = true;
     break;
    }
   } // loop on the c ops of the configuration

   // if previous loop is broken, partial_trace is 0 
   auto partial_trace_no_emin = (loop_broken ? 0 : trace(U_boundary.transpose() * M));
   //auto partial_trace_no_emin = (loop_broken ? 0 : trace(U_boundary.transpose() * M_work(range(0,space_dim),range(0,space_dim))));
   auto partial_trace = partial_trace_no_emin * exp_no_emin;
   if (std::abs(partial_trace_no_emin) > space_dim * 1.0000001) throw "halte la !"; // CHECK conjecture

   //std::cout  << bl << " "<< partial_trace << " "<<full_trace << " bound = "<< std::exp(-to_sort[bl].first) << std::endl;
   if (bl == 0) first_term = partial_trace;
   full_trace += partial_trace;
   partial_trace_of_block[block_index] = partial_trace;
   for (int j = 0; j < trunc_block.size(); ++j){
    if (bl < trunc_block[j]) partial_trace_up_to_block[j] += partial_trace;
   }

  } // -.-.-.-.-.-.-.-.-.  choice of trace computation method -.-.-.-.-.-.-

  // DEBUG
  // if (bl ==5) break;
 } // end of loop on blocks
 // std::cout  << "-----------"<< std::endl;

 if (make_histograms) {
  auto abs_trace = std::abs(full_trace);
  if (abs_trace > 0) histos["FirsTerm_FullTrace"] << std::abs(first_term) / abs_trace;
  histos["FullTrace_ExpSumMin"] << std::abs(full_trace) / std::exp(-to_sort[0].first);
  for (int i = 0; i < trunc_block.size() ; ++i){
   std::stringstream s;
   s << "TruncatedTrace_over_FullTrace" << trunc_block[i];
   if (abs_trace > 0) histos[s.str()] << std::abs(partial_trace_up_to_block[i]) / abs_trace;
  }
  histo_bs_block << to_sort[0].second;
 }
 for (int i = 0; i < partial_trace_of_block.size(); ++i) partial_over_full_trace[i] += partial_trace_of_block[i] / full_trace;
 return full_trace;
}

// ------------------- EXPERIMENTAL full trace computation -------------
// try a direct computation, no estimator, pushing always the largest partial path... 
struct w_pt { 
 int time_index;
 long block_index, block_start;
};

atomic_correlators_worker::result_t atomic_correlators_worker::full_trace2(double epsilon) {

 int n_blocks = sosp.n_subspaces();
 int config_size = config->size();
 std::vector<result_t> partial_trace_of_block(n_blocks,0);
 std::vector<result_t> partial_trace_up_to_block(trunc_block.size(),0);

 auto config_table = make_config_table(config);

 //std::cout  << * config << std::endl; 
 if (config_size == 0) { // special case in empty configuration : all paths are already finished !
  std::cout  << "empty config " << std::endl;
  double res = 0;
  for (int n = 0; n < n_blocks; ++n) {
   auto d = sosp.subspace(n).dimension();
   for (int k = 0; k < d; ++k)
    res += std::exp(- config->beta() * (sosp.get_eigensystems()[n].eigenvalues(k))); // - eigensystem.eigenvalues(0)));
  }
  return res;
 }

 auto dtau0 = double(config->lowest_time_operator()->first);
 double w_sum = 0;
 result_t full_trace = 0;

 //int max_subspace_dim=0;
 //for (int nsp = 0; nsp < sosp.n_subspaces(); ++nsp) max_subspace_dim = std::max(max_subspace_dim, sosp.subspace(nsp).dimension());
 
 auto m_work = std::vector<triqs::arrays::matrix<double>>(n_blocks); //, triqs::arrays::matrix<double>{max_subspace_dim,max_subspace_dim});
 std::multimap<double, w_pt, std::greater<double>> w_p_path; // weighted partial paths

 // init the paths
 for (int n = 0; n < n_blocks; ++n) {
  auto d = sosp.subspace(n).dimension();
  m_work[n] = sosp.get_eigensystems()[n].unitary_matrix;
  for (int k = 0; k < d; ++k)
   m_work[n](k, range()) *= std::exp(-dtau0 * (sosp.get_eigensystems()[n].eigenvalues(k))); // - eigensystem.eigenvalues(0)));
  double w = norm1(m_work[n]);
  w_sum += w;
  w_p_path.insert({w, {0, n, n}});
 }

 // main loop : continue the path with the largest weight
 while ((w_p_path.size()) && (w_sum > epsilon * std::abs(full_trace))) {

  // pick up the the path with the largest weight
  auto p = w_p_path.begin();
  w_p_path.erase(p); // always erase the point what we will now treat
  auto w = p->first;
  w_sum -= w;
  auto B_start = p->second.block_start;
  auto Bi = p->second.block_index;
  auto i = p->second.time_index;

  // act with c op
  auto Bf = sosp.fundamental_operator_connect_from_linear_index(config_table[i].dag, config_table[i].n, Bi);
  if (Bf == -1) continue; // the path has ended ...

  // apply operator
  auto const& Cop = sosp.fundamental_operator_matrix_from_linear_index(config_table[i].dag, config_table[i].n, Bi);
  auto const& eigensystem_i = sosp.get_eigensystems()[Bi];
  auto const& eigensystem_f = sosp.get_eigensystems()[Bf];
  auto space_dim_i = eigensystem_i.eigenvalues.size();
  auto space_dim_f = eigensystem_f.eigenvalues.size();

  if ((space_dim_i == 1) && (space_dim_f == 1)) // some quick optimisation
   m_work[B_start] *= Cop(0, 0);
  else {
   auto M2 = Cop * m_work[B_start];
   swap(m_work[B_start], M2);
  }

  // apply exponential.
  for (int n = 0; n < space_dim_f; ++n)
   m_work[B_start](n, range()) *= std::exp(-config_table[i].dtau * (eigensystem_f.eigenvalues(n))); //- eigensystem.eigenvalues(0)));

  // recompute the weight
  auto w2 = norm1(m_work[B_start]);

  if (i+1 == config_size) { // the path is at the end
   auto const& U_boundary = sosp.get_eigensystems()[B_start].unitary_matrix;
   auto partial_trace = trace(U_boundary.transpose() * m_work[B_start]);
   full_trace += partial_trace;
  } else { // the path is not at the end
   w_p_path.insert({w2, {i+1, Bf, B_start}});
   w_sum += w2;
  }

 } // end main loop

 // analysis
 if (make_histograms) {
 }

 return full_trace;
}

}// namespace

