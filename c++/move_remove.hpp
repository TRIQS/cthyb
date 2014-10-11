/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once
#include <algorithm>
#include "./qmc_data.hpp"

namespace cthyb {

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
  std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  std::cerr << "* Attempt for move_remove_c_cdag (block " << block_index << ")" << std::endl;
  std::cerr << "* Configuration before:" << std::endl << config;
  data.imp_trace.tree.graphviz(std::ofstream("tree_before"));
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
  tau1 = data.imp_trace.try_delete(num_c, block_index, false);
  tau2 = data.imp_trace.try_delete(num_c_dag, block_index, true);

  auto det_ratio = det.try_remove(num_c_dag, num_c);

  // acceptance probability
  double t_ratio = std::pow(block_size * config.beta() / double(det_size), 2);
  
  // For quick abandon 
  double random_number = rng.preview();
  if (random_number == 0.0) return 0;
  double p_yee = std::abs(det_ratio / t_ratio/ data.trace);

  // recompute the trace
  new_trace = data.imp_trace.estimate(p_yee, random_number);
  if (new_trace == 0.0) return 0;
  auto trace_ratio = new_trace / data.trace;
  if (!std::isfinite(trace_ratio)) TRIQS_RUNTIME_ERROR << "trace_ratio not finite" << new_trace << "  "<< data.trace<<"  "<< new_trace /data.trace ;
 
  mc_weight_type p = trace_ratio * det_ratio;

#ifdef EXT_DEBUG
  std::cerr << "Trace ratio: " << trace_ratio << '\t';
  std::cerr << "Det ratio: " << det_ratio << '\t';
  std::cerr << "Prefactor: " << t_ratio << '\t';
  std::cerr << "Weight: " << p/ t_ratio << std::endl;
#endif

  if (!std::isfinite(p)) TRIQS_RUNTIME_ERROR << "(remove) p not finite :" << p;
  if (!std::isfinite(p / t_ratio)) TRIQS_RUNTIME_ERROR << "p / t_ratio not finite p : " << p << " t_ratio :  "<< t_ratio;
  return p / t_ratio;
 }

 //----------------

 mc_weight_type accept() {

  // remove from the configuration
  config.erase(tau1);
  config.erase(tau2);

  // remove in the cache tree  
  data.imp_trace.confirm_delete();
  
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
  data.imp_trace.cancel_delete();
#ifdef EXT_DEBUG
  std::cerr << "* Configuration after: " << std::endl;
  std::cerr << config;
#endif

 }
};
}
