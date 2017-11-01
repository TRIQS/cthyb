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

#include "./global.hpp"

namespace cthyb {

  move_global::move_global(std::string const &name, indices_map_t const &substitution_map, qmc_data &data, mc_tools::random_generator &rng)
     : name(name),
       data(data),
       config(data.config),
       rng(rng),
       substitute_c(data.linindex.size()),
       substitute_c_dag(data.linindex.size()),
       x(data.dets.size()),
       y(data.dets.size()) {

    auto const &fops = data.h_diag.get_fops();

    // Inverse of data.linindex
    std::vector<std::pair<int, int>> lin_to_block_inner(data.linindex.size());
    for (auto const &l : data.linindex) lin_to_block_inner[l.second] = l.first;

    bool identity = true;
    for (int lin = 0; lin < lin_to_block_inner.size(); ++lin) {
      int new_lin, new_block, new_inner;

      // Does operator with linear index lin have a mapping in substitution_map?
      auto it = std::find_if(std::begin(substitution_map), std::end(substitution_map),
                             [&fops, lin](indices_map_t::value_type const &kv) { return fops[kv.first] == lin; });

      // If it does not, it is substituted by itself (subst_linear = lin)
      new_lin = (it != std::end(substitution_map)) ? fops[it->second] : lin;
      std::tie(new_block, new_inner) = lin_to_block_inner[new_lin];

      if (new_lin != lin) {
        identity = false;
        affected_blocks.insert(lin_to_block_inner[lin].first);
        affected_blocks.insert(new_block);
      }

      substitute_c[lin]     = op_desc{new_block, new_inner, false, new_lin};
      substitute_c_dag[lin] = op_desc{new_block, new_inner, true, new_lin};
    }

    if (identity) std::cerr << "WARNING: global move '" << name << "' changes no operator indices, therefore is useless." << std::endl;
  }

  mc_weight_t move_global::attempt() {

#ifdef EXT_DEBUG
    std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cerr << "In config " << config.get_id() << std::endl;
    std::cerr << "* Attempt for move move_global (" << name << ")" << std::endl;
#endif

    updated_ops.clear();
    for (auto const &o : data.config) {
      auto const &tau    = o.first;
      auto const &old_op = o.second;
      auto const &new_op = (old_op.dagger ? substitute_c_dag : substitute_c)[old_op.linear_index];
      if (old_op.linear_index != new_op.linear_index) updated_ops.emplace(tau, new_op);
    }

#ifdef EXT_DEBUG
    std::cerr << updated_ops.size() << " out of " << data.config.size() << " operators can be changed" << std::endl;
#endif

    // No operators can be updated...
    if (!updated_ops.size()) return 0;

    // Choose a random number of operators, which will not actually be updated
    // (we always update at least one operator)
    int n_no_update = rng(updated_ops.size());
    // Remove some operators
    for (int i = 0; i < n_no_update; ++i) {
      auto it = std::begin(updated_ops);
      std::advance(it, rng(updated_ops.size()));
      updated_ops.erase(it);
    }

#ifdef EXT_DEBUG
    std::cerr << updated_ops.size() << " operators will actually be changed" << std::endl;
#endif

    // Derive new arguments of the dets
    for (auto block_index : affected_blocks) {
      x[block_index].clear();
      y[block_index].clear();
    }

    for (auto const &o : data.config) {
      auto const &tau    = o.first;
      auto it            = updated_ops.find(tau);
      auto const &new_op = it == updated_ops.end() ? o.second : it->second;
      (new_op.dagger ? x : y)[new_op.block_index].emplace_back(tau, new_op.inner_index);
    }

    for (auto block_index : affected_blocks) {
      // New configuration is not compatible with gf_struct
      if (x[block_index].size() != y[block_index].size()) return 0;
    }

    // Try refill determinants
    mc_weight_t det_ratio = 1;
    for (auto block_index : affected_blocks) {
      auto &det                   = data.dets[block_index];
      mc_weight_t block_det_ratio = det.try_refill(x[block_index], y[block_index]);
      if (block_det_ratio == .0) {
#ifdef EXT_DEBUG
        std::cerr << "block_det_ratio[" << block_index << "] = 0" << std::endl;
#endif
        return 0;
      }
      det_ratio *= block_det_ratio;
    }

    // For quick abandon
    double random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double p_yee = std::abs(det_ratio / data.atomic_weight);

    data.imp_trace.try_replace(updated_ops);

    // computation of the new trace after insertion
    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) {
#ifdef EXT_DEBUG
      std::cerr << "trace == 0" << std::endl;
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
    std::cerr << "p_yee: " << p_yee << std::endl;
    std::cerr << "Weight: " << p << std::endl;
#endif

    return p;
  }

  mc_weight_t move_global::accept() {

    for (auto const &o : updated_ops) data.config.replace(o.first, o.second);
    config.finalize();

    for (auto block_index : affected_blocks) data.dets[block_index].complete_operation();

    data.update_sign();
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;

    data.imp_trace.confirm_replace();

#ifdef EXT_DEBUG
    std::cerr << "* Move move_global '" << name << "' accepted" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    for (int block_index : affected_blocks) check_det_sequence(data.dets[block_index], config.get_id());
#endif

    return data.current_sign / data.old_sign;
  }

  void move_global::reject() {

    config.finalize();
    data.imp_trace.cancel_replace();
    for (auto block_index : affected_blocks) data.dets[block_index].reject_last_try();

#ifdef EXT_DEBUG
    std::cerr << "* Move move_global '" << name << "' rejected" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    for (int block_index : affected_blocks) check_det_sequence(data.dets[block_index], config.get_id());
#endif
  }
}
