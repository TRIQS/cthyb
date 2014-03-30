#pragma once
#include <algorithm>
#include "./qmc_data.hpp"

namespace cthyb_matrix {

// Removal of C, C^dagger operator
class move_remove_c_cdag {

 qmc_data& data;
 configuration& config;
 mc_tools::random_generator& rng;
 int block_index, block_size;
 qmc_data::trace_t new_trace;
 time_pt tau1, tau2;

 public:
 
 //----------------------------------

 move_remove_c_cdag(int block_index, int block_size, qmc_data& data, mc_tools::random_generator& rng)
    : data(data), config(data.config), rng(rng), block_index(block_index), block_size(block_size) {}

 //----------------
 
 mc_weight_type attempt() {

#ifdef EXT_DEBUG
  config.id++;
  config.print_debug();
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "* Attempt for move_remove_c_cdag (block " << block_index << ")" << std::endl;
  std::cerr << "* Configuration before:" << std::endl << config;
  data.atomic_corr.tree.graphviz(std::ofstream("tree_before"));
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

  // now mark 2 nodes for deletion
  tau1 = data.atomic_corr.soft_delete_n_th_operator(num_c, block_index, false);
  tau2 = data.atomic_corr.soft_delete_n_th_operator(num_c_dag, block_index, true);

  new_trace = data.atomic_corr.estimate();
  if (new_trace == 0.0) return 0;
  auto trace_ratio = new_trace / data.trace;

  if (!std::isfinite(trace_ratio)) TRIQS_RUNTIME_ERROR << "trace_ratio not finite" << new_trace << "  "<< data.trace<<"  "<< new_trace /data.trace ;
  if (!std::isnormal(trace_ratio)) TRIQS_RUNTIME_ERROR << "trace_ratio not normal" << new_trace << "  "<< data.trace <<"  "<< new_trace /data.trace ;
 
  auto det_ratio = det.try_remove(num_c_dag, num_c);

  // acceptance probability
  mc_weight_type p = trace_ratio * det_ratio;
  double t_ratio = std::pow(block_size * config.beta() / double(det_size), 2);

#ifdef EXT_DEBUG
  std::cerr << "Trace ratio: " << trace_ratio << '\t';
  std::cerr << "Det ratio: " << det_ratio << '\t';
  std::cerr << "Prefactor: " << t_ratio << '\t';
  std::cerr << "Weight: " << p* t_ratio << std::endl;
#endif

  if (!std::isfinite(p)) TRIQS_RUNTIME_ERROR << "(remove) p not finite :" << p;
  if (!std::isfinite(p * t_ratio)) TRIQS_RUNTIME_ERROR << "p * t_ratio not finite" << p* t_ratio;
  return p / t_ratio;
 }

 //----------------

 mc_weight_type accept() {

  // remove from the configuration
  config.erase(tau1);
  config.erase(tau2);

  // remove in the cache tree  
  data.atomic_corr.confirm_soft_delete();
  
  // remove in the config
  data.dets[block_index].complete_operation();
  data.update_sign();
  data.trace = new_trace;
#ifdef EXT_DEBUG
  std::cerr << "* Configuration after: " << std::endl;
  std::cerr << config;
#endif

  return data.current_sign / data.old_sign;
 }

 //----------------

 void reject() {
  data.atomic_corr.clean_soft_delete();
#ifdef EXT_DEBUG
  std::cerr << "* Configuration after: " << std::endl;
  std::cerr << config;
#endif

 }
};
}
