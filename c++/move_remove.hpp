#pragma once
#include <algorithm>
#include "./qmc_data.hpp"

namespace cthyb_krylov {

using mc_weight_type = std::complex<double>;

// Removal of C, C^dagger operator
class move_remove_c_cdag {

 qmc_data& data;
 configuration& config;
 mc_tools::random_generator& rng;

 int block_index, block_size;
 qmc_data::trace_t new_trace;

 std::vector<std::pair<configuration::oplist_t::key_type, configuration::oplist_t::mapped_type>> removed_ops;

 // remove the n-th operator of type dag. Returns the time_pt where it was removed
 time_pt remove_op(int n, bool dag) {
  // Find the position of the operator to remove
  int i = 0;
  auto it = std::find_if(config.begin(), config.end(), [&](typename configuration::oplist_t::value_type const& op) {
   if (op.second.dagger == dag && op.second.block_index == block_index) ++i;
   return i == n + 1;
  });
  assert(it != config.end());
  removed_ops.push_back(*it); // store the pair (time, operator number) which was here, to put it back in reject if necessary
  auto t = it->first;
  config.erase(it);
  return t;
 }

 public:
 //----------------------------------

 move_remove_c_cdag(int block_index, int block_size, qmc_data& data, mc_tools::random_generator& rng)
    : data(data), config(data.config), rng(rng), block_index(block_index), block_size(block_size) {}

 //----------------
 mc_weight_type attempt() {

#ifdef EXT_DEBUG
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "* Attempt for move_remove_c_cdag (block " << block_index << ")" << std::endl;
  std::cerr << "* Configuration before:" << std::endl << config;
#endif

  // the det has to be recomputed each time, since global moves will change it
  auto& det = data.dets[block_index];

  // Pick up a couple of C, Cdagger to remove at random
  // Remove the operators from the traces
  int det_size = det.size();
  if (det_size == 0) return 0; // nothing to remove
  int num_c_dag = rng(det_size), num_c = rng(det_size);
#ifdef EXT_DEBUG
  std::cerr << "* Proposing to remove: ";
  std::cerr << num_c_dag << "-th Cdag(" << block_index << ",...), ";
  std::cerr << num_c << "-th C(" << block_index << ",...)" << std::endl;
#endif
  auto t1 = remove_op(num_c, false);
  auto t2 = remove_op(num_c_dag, true);

  new_trace = data.atomic_corr.estimate(t1, t2);
  if (new_trace == 0.0) return 0;
  auto trace_ratio = new_trace / data.trace;

  auto det_ratio = det.try_remove(num_c_dag, num_c);

  // acceptance probability
  mc_weight_type p = trace_ratio * det_ratio;
  double t_ratio = std::pow(block_size * config.beta() / double(det_size), 2);

#ifdef EXT_DEBUG
  std::cerr << "Trace ratio: " << trace_ratio << '\t';
  std::cerr << "Det ratio: " << det_ratio << '\t';
  std::cerr << "Prefactor: " << t_ratio << '\t';
  std::cerr << "Weight: " << p* t_ratio << std::endl;
  std::cerr << "* Configuration after: " << std::endl;
  std::cerr << config;
#endif

  return p / t_ratio;
 }

 //----------------

 mc_weight_type accept() {
  removed_ops.clear();
  data.dets[block_index].complete_operation();
  data.update_sign();
  data.trace = new_trace;
  data.atomic_corr.cache_update();
  //for (int i =0; i<100; ++i) data.atomic_corr.cache_update();
  return data.current_sign / data.old_sign;
 }

 //----------------

 void reject() {
  for (auto& p : removed_ops) config.insert(p);
  removed_ops.clear();
 }
};
}
