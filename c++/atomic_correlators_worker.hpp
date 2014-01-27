#pragma once
#include "configuration.hpp"
#include "sorted_spaces.hpp"
#include "exp_h_worker.hpp"
#include <triqs/parameters.hpp>

#include "triqs/statistics/histograms.hpp"
namespace cthyb_krylov {

/**
 * A worker that computes the trace using krylov method, for a given configuration.
 * Has to live longer than the configuration...
 */
class atomic_correlators_worker {
 public:
 using result_t = std::complex<double>;

 atomic_correlators_worker(configuration& c, sorted_spaces const& sosp_, utility::parameters const &p);

 ~atomic_correlators_worker();

 sorted_spaces const& get_sorted_spaces() const { return sosp; }

 result_t estimate(time_pt tau1, time_pt tau2);
 result_t estimate();
 result_t full_trace_over_estimator();

 void cache_update();

 private:

 // Various possible algorithms
 enum class estimator_method_t {None, Simple, WithCache};

 struct cache_point_t {
  struct pt_t {
   int current_block_number; // -1 means no block
   double emin_dtau_acc;
  };
  std::vector<pt_t> r, l; // precomputation at the right and the left of the point
  cache_point_t(int n_blocks) : r(n_blocks, {-1, 0}), l(n_blocks, {-1, 0}) {}
 };

 cache_point_t make_cache_point() const { return {sosp.n_subspaces()};}

 using cache_t = std::map<time_pt, cache_point_t, std::greater<time_pt>>;
 // using oplist_t=boost::container::flat_map<time_pt, cache_point_t, std::greater<time_pt>> ;
 cache_t cache;

 const configuration* config;                                     // must exists longer than this object.
 sorted_spaces const& sosp;                                       // The sorted space

 using state_t = state<sub_hilbert_space, double, false>;
 exp_h_worker<imperative_operator<sub_hilbert_space, false>, state_t> exp_h;

 bool make_histograms;                                            // Do we make the Histograms ?
 std::map<std::string, statistics::histogram_segment_bin> histos; // Analysis histograms
 statistics::histogram histo_bs_block;                            // Histogram of the boundary state
 // statistics::histogram histo_opcount;                             // Histogram of number of operators in non-zero path
 statistics::histogram histo_trace_null_struc;

 std::vector<statistics::histogram> histo_n_blocks_after_steps;
 std::vector<statistics::histogram> histo_opcount;
 estimator_method_t estimator_method;
 bool use_truncation;
 bool use_old_trace;
 std::vector<result_t> time_spent_in_block;
 std::vector<result_t> partial_over_full_trace;

 result_t estimate_with_cache(time_pt tau1, time_pt tau2);
 result_t estimate_simple(bool no_exp = false);
 result_t full_trace();
};
}
