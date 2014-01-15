#pragma once
#include "./atomic_correlators_worker.hpp"
#include <triqs/gfs.hpp>
#include <triqs/det_manip.hpp>
#include <triqs/utility/serialization.hpp>

namespace cthyb_krylov {
using namespace triqs::gfs;

/**
 * The data of the Monte carlo
 */
struct qmc_data {

 configuration config;                  // Configuration
 sorted_spaces const &sosp;             // Diagonalization of the atomic problem
 mutable atomic_correlators_worker atomic_corr; // Calculator of the trace

 using trace_t = atomic_correlators_worker::result_t;

 /// This callable object adapts the Delta function for the call of the det.
 struct delta_block_adaptor {
  gf_const_view<imtime> delta_block;

  delta_block_adaptor(gf_const_view<imtime> const &delta_block) : delta_block(delta_block) {}
  delta_block_adaptor(delta_block_adaptor const &) = default;
  delta_block_adaptor(delta_block_adaptor &&) = default;
  delta_block_adaptor &operator=(delta_block_adaptor const &) = delete; // forbid assignment
  delta_block_adaptor &operator=(delta_block_adaptor &&a) = default;

  // no need of argument_type, return_type : det_manip now synthetize everything (need to UPDATE doc).
  double operator()(std::pair<time_pt, int> const &x, std::pair<time_pt, int> const &y) const {
   double res = delta_block[closest_mesh_pt(double(x.first - y.first))](x.second, y.second);
   return (x.first >= y.first ? res : -res); // x,y first are time_pt, wrapping is automatic in the - operation, but need to
                                             // compute the sign
  }
 };

 std::vector<det_manip::det_manip<delta_block_adaptor>> dets; // The determinants
 int current_sign, old_sign;                                  // Permutation prefactor
 trace_t trace;                                               // The current value of the trace
 bool use_quick_trace_estimator;

 // construction and the basics stuff. value semantics, except = ?
 qmc_data(utility::parameters const &p, sorted_spaces const &sosp, block_gf_const_view<imtime> delta)
    : config(p["beta"]),
      sosp(sosp),
      atomic_corr(config, sosp, p["krylov_gs_energy_convergence"], p["krylov_small_matrix_size"], p["make_path_histograms"],
                  p["use_quick_trace_estimator"], p["trace_estimator_n_blocks_guess"], p["use_truncation"], p["use_old_trace"]),
      current_sign(1),
      old_sign(1) {
  use_quick_trace_estimator = p["use_quick_trace_estimator"];
  trace = atomic_corr();
  dets.clear();
  for (auto const &bl : delta.mesh()) dets.emplace_back(delta_block_adaptor(delta[bl]), 100);
 }

 qmc_data(qmc_data const &) = default;
 qmc_data &operator=(qmc_data const &) = delete;

 void update_sign() {

  int s = 0;
  size_t num_blocks = dets.size();
  vector<int> n_op_with_a_equal_to(num_blocks, 0), n_ndag_op_with_a_equal_to(num_blocks, 0);

  // In this first part we compute the sign to bring the configuration to
  // d^_1 d^_1 d^_1 ... d_1 d_1 d_1   d^_2 d^_2 ... d_2 d_2   ...   d^_n .. d_n

  // loop over the operators "op" in the trace (right to left)
  for (auto const &op : config.oplist) {

   // how many operators with an 'a' larger than "op" are there on the left of "op"?
   for (int a = op.second.block_index + 1; a < num_blocks; ++a) s += n_op_with_a_equal_to[a];
   n_op_with_a_equal_to[op.second.block_index]++;

   // if "op" is not a dagger how many operators of the same a but with a dagger are there on his right?
   if (op.second.dagger)
    s += n_ndag_op_with_a_equal_to[op.second.block_index];
   else
    n_ndag_op_with_a_equal_to[op.second.block_index]++;
  }

  // Now we compute the sign to bring the configuration to
  // d_1 d^_1 d_1 d^_1 ... d_1 d^_1   ...   d_n d^_n ... d_n d^_n
  for (int block_index = 0; block_index < num_blocks; block_index++) {
   int n = dets[block_index].size();
   s += n * (n + 1) / 2;
  }

  old_sign = current_sign;
  current_sign = (s % 2 == 0 ? 1 : -1);
 }

 template <class Archive> void serialize(Archive &ar, const unsigned int version) {
  ar &boost::serialization::make_nvp("configuration", config) & boost::serialization::make_nvp("dets", dets) &
      boost::serialization::make_nvp("atomic_corr", atomic_corr) & boost::serialization::make_nvp("old_sign", old_sign) &
      boost::serialization::make_nvp("current_sign", current_sign);
 }
};
}

