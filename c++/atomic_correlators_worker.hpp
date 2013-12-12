#pragma once
#include "configuration.hpp"
#include "sorted_spaces.hpp"
#include "exp_h_worker.hpp"

#include "./histograms.hpp"
namespace cthyb_krylov {

/**
 * A worker that computes the trace using krylov method, for a given configuration.
 * Has to live longer than the configuration...
 */
class atomic_correlators_worker {
 public:
 
 using result_t = std::complex<double>;
 
 atomic_correlators_worker(configuration& c, sorted_spaces const& sosp_, double gs_energy_convergence, int small_matrix_size);

 result_t operator()(); // recompute and return the full trace

 const sorted_spaces& get_sorted_spaces() const { return sosp; }
 
 private:
 const configuration* config; // must exists longer than this object.
 sorted_spaces sosp;          // The sorted space
 int small_matrix_size;// The minimal size of a matrix to be treated with exp_h_matrix
 std::map<std::string,statistics::histogram_segment_bin> histos;
 statistics::histogram histo_bs_block; 

 using state_t = state<sub_hilbert_space, double, false>;
 exp_h_worker<imperative_operator<sub_hilbert_space, false>, state_t> exp_h;
};
}
