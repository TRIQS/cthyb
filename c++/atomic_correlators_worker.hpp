#pragma once
#include "configuration.hpp"
#include "sorted_spaces.hpp"
#include "exp_h_worker.hpp"

#include "triqs/statistics/histograms.hpp"
namespace cthyb_krylov {

/**
 * A worker that computes the trace using krylov method, for a given configuration.
 * Has to live longer than the configuration...
 */
class atomic_correlators_worker {
 public:
 using result_t = std::complex<double>;

 atomic_correlators_worker(configuration& c, sorted_spaces const& sosp_, double gs_energy_convergence, int small_matrix_size,
                           bool make_histograms, bool use_quick_trace_estimator, int trace_estimator_n_blocks_guess,
                           bool use_truncation, bool use_old_trace);

 ~atomic_correlators_worker();

 result_t operator()(); // recompute and return the full trace
 result_t full_trace_over_estimator();

 sorted_spaces const& get_sorted_spaces() const { return sosp; }

 private:
 
 struct cache_point_t {
  struct pt_t {
   int current_block_number;//, start_block_number;
   double emin_dtau_acc;
  };
  std::vector<pt_t> r, l; // precomputation at the right and the left of the point
  cache_point_t(int n_blocks) : r{n_blocks, {0, 0}}, l{n_blocks, {0, 0}} {}
  //cache_point_t(int n_blocks) : r{n_blocks, {0, 0, 0}}, l{n_blocks, {0, 0, 0}} {}
 };

 using cache_t = std::map<time_pt, cache_point_t, std::greater<time_pt>>;
 // using oplist_t=boost::container::flat_map<time_pt, cache_point_t, std::greater<time_pt>> ;
 cache_t cache;

 const configuration* config;                                     // must exists longer than this object.
 sorted_spaces const& sosp;                                       // The sorted space
 int small_matrix_size;                                           // The minimal size of a matrix to be treated with exp_h_matrix
 bool make_histograms;                                            // Do we make the Histograms ?
 std::map<std::string, statistics::histogram_segment_bin> histos; // Analysis histograms
 statistics::histogram histo_bs_block;                            // Histogram of the boundary state
 //statistics::histogram histo_opcount;                             // Histogram of number of operators in non-zero path
 statistics::histogram histo_trace_null_struc; 

 std::vector<statistics::histogram> histo_n_blocks_after_steps;
 std::vector<statistics::histogram> histo_opcount;
 bool use_quick_trace_estimator;
 int trace_estimator_n_blocks_guess;
 bool use_truncation;
 bool use_old_trace;
 std::vector<result_t> time_spent_in_block;
 std::vector<result_t> partial_over_full_trace;
 triqs::arrays::array<int,2> block_died_anal;
 triqs::arrays::array<int,2> block_died_num;

 using state_t = state<sub_hilbert_space, double, false>;
 exp_h_worker<imperative_operator<sub_hilbert_space, false>, state_t> exp_h;

 result_t estimate();
 result_t full_trace();
};
}
