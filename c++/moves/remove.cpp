/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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

#include "./remove.hpp"

namespace cthyb {

// -----------------------------------------------------------------------------
histogram* move_remove_c_cdag::add_histo(std::string const& name, histo_map_t* histos) {
 if (!histos) return nullptr;
 auto new_histo = histos->insert({name, {.0, config.beta(), 100}});
 return &(new_histo.first->second);
}

// -----------------------------------------------------------------------------
mc_weight_t move_remove_c_cdag::attempt() {

#ifdef EXT_DEBUG
 std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
 std::cerr << "In config " << config.get_id() << std::endl;
 std::cerr << "* Attempt for move_remove_c_cdag (block " << block_index << ")" << std::endl;
#endif

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

 // record the length of the proposed removal
 dtau = double(tau2 - tau1);
 if (histo_proposed) *histo_proposed << dtau;

 auto det_ratio = det.try_remove(num_c_dag, num_c);

 // proposition probability
 auto t_ratio = std::pow(block_size * config.beta() / double(det_size), 2); // Size of the det before the try_delete!

 // For quick abandon
 double random_number = rng.preview();
 if (random_number == 0.0) return 0;
 double p_yee = std::abs(det_ratio / t_ratio / data.atomic_weight);

 // recompute the atomic_weight
 std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
 if (new_atomic_weight == 0.0) {
#ifdef EXT_DEBUG
  std::cerr << "atomic_weight == 0" << std::endl;
#endif
  return 0;
 }
 auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
 if (!isfinite(atomic_weight_ratio))
  TRIQS_RUNTIME_ERROR << "atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight << " "
                      << new_atomic_weight / data.atomic_weight << " in config " << config.get_id();

 mc_weight_t p = atomic_weight_ratio * det_ratio;

#ifdef EXT_DEBUG
 std::cerr << "Trace ratio: " << atomic_weight_ratio << '\t';
 std::cerr << "Det ratio: " << det_ratio << '\t';
 std::cerr << "Prefactor: " << t_ratio << '\t';
 std::cerr << "Weight: " << p / t_ratio << std::endl;
#endif

 if (!isfinite(p)) TRIQS_RUNTIME_ERROR << "(remove) p not finite :" << p << " in config " << config.get_id();
 if (!isfinite(p / t_ratio))
  TRIQS_RUNTIME_ERROR << "p / t_ratio not finite p : " << p << " t_ratio :  " << t_ratio << " in config " << config.get_id();
 return p / t_ratio;
}

// -----------------------------------------------------------------------------
mc_weight_t move_remove_c_cdag::accept() {

 // remove from the tree
 data.imp_trace.confirm_delete();

 // remove from the configuration
 config.erase(tau1);
 config.erase(tau2);
 config.finalize();

 // remove from the determinants
 data.dets[block_index].complete_operation();
 data.update_sign();
 data.atomic_weight = new_atomic_weight;
 data.atomic_reweighting = new_atomic_reweighting;
 if (histo_accepted) *histo_accepted << dtau;

#ifdef EXT_DEBUG
 std::cerr << "* Move move_remove_c_cdag accepted" << std::endl;
 std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
 check_det_sequence(data.dets[block_index], config.get_id());
#endif

 return data.current_sign / data.old_sign;
}

// -----------------------------------------------------------------------------
void move_remove_c_cdag::reject() {

 config.finalize();
 data.imp_trace.cancel_delete();

#ifdef EXT_DEBUG
 std::cerr << "* Move move_remove_c_cdag rejected" << std::endl;
 std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
 check_det_sequence(data.dets[block_index], config.get_id());
#endif
}
}
