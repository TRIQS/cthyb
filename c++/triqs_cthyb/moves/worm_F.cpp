/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2026
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

#include "./worm_F.hpp"

#include <cmath>

namespace triqs_cthyb {

  namespace {
    double worm_phase_space(qmc_data const &data, double worm_eta) {
      return worm_eta * std::pow(data.config.beta(), 2) * double(data.worm.n_components());
    }
  } // namespace

  move_worm_insert_F::move_worm_insert_F(qmc_data &data, mc_tools::random_generator &rng, double worm_eta)
     : data(data), rng(rng), worm_eta(worm_eta) {}

  mc_weight_t move_worm_insert_F::attempt() {
    if (data.worm.in_F() || data.worm.components.empty()) return 0;

    component_index      = rng(data.worm.n_components());
    auto const &component = data.worm.components[component_index];
    tau_Q                = data.tau_seg.get_random_pt(rng);
    tau_cdag             = data.tau_seg.get_random_pt(rng);

    try {
      data.imp_trace.try_insert(tau_Q, *data.worm.Q_ops[component.block][component.inner_Q]);
      data.imp_trace.try_insert(tau_cdag, data.worm.cdag_ops[component.block][component.inner_cdag]);
    } catch (rbt_insert_error const &) {
      data.imp_trace.cancel_insert();
      return 0;
    }

    auto t_ratio = worm_phase_space(data, worm_eta);

    double random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double p_yee = std::abs(t_ratio / data.atomic_weight);

    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) return 0;

    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(worm insert F) atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight;

    return atomic_weight_ratio * t_ratio;
  }

  mc_weight_t move_worm_insert_F::accept() {
    auto const &component = data.worm.components[component_index];

    data.imp_trace.confirm_insert();
    data.worm.sector          = qmc_data::worm_data_t::sector_t::F;
    data.worm.component_index = component_index;
    data.worm.block           = component.block;
    data.worm.inner_Q         = component.inner_Q;
    data.worm.inner_cdag      = component.inner_cdag;
    data.worm.tau_Q           = tau_Q;
    data.worm.tau_cdag        = tau_cdag;
    data.atomic_weight        = new_atomic_weight;
    data.atomic_reweighting   = new_atomic_reweighting;
    return 1;
  }

  void move_worm_insert_F::reject() { data.imp_trace.cancel_insert(); }

  move_worm_remove_F::move_worm_remove_F(qmc_data &data, mc_tools::random_generator &rng, double worm_eta) : data(data), rng(rng), worm_eta(worm_eta) {}

  mc_weight_t move_worm_remove_F::attempt() {
    if (data.worm.in_Z()) return 0;

    data.imp_trace.try_delete(data.worm.tau_Q);
    data.imp_trace.try_delete(data.worm.tau_cdag);

    auto t_ratio = worm_phase_space(data, worm_eta);
    double random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double p_yee = std::abs(1.0 / (t_ratio * data.atomic_weight));

    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) return 0;

    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(worm remove F) atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight;

    return atomic_weight_ratio / t_ratio;
  }

  mc_weight_t move_worm_remove_F::accept() {
    data.imp_trace.confirm_delete();
    data.worm.sector          = qmc_data::worm_data_t::sector_t::Z;
    data.worm.component_index = -1;
    data.worm.block           = -1;
    data.worm.inner_Q         = -1;
    data.worm.inner_cdag      = -1;
    data.atomic_weight        = new_atomic_weight;
    data.atomic_reweighting   = new_atomic_reweighting;
    return 1;
  }

  void move_worm_remove_F::reject() { data.imp_trace.cancel_delete(); }

  move_worm_shift_F::move_worm_shift_F(qmc_data &data, mc_tools::random_generator &rng) : data(data), rng(rng) {}

  mc_weight_t move_worm_shift_F::attempt() {
    if (data.worm.in_Z()) return 0;

    shift_Q = (rng(2) == 0);
    tau_old = shift_Q ? data.worm.tau_Q : data.worm.tau_cdag;
    tau_new = data.tau_seg.get_random_pt(rng);
    op_new  = shift_Q ? *data.worm.Q_ops[data.worm.block][data.worm.inner_Q] : data.worm.cdag_ops[data.worm.block][data.worm.inner_cdag];

    data.imp_trace.try_delete(tau_old);
    try {
      data.imp_trace.try_insert(tau_new, op_new);
    } catch (rbt_insert_error const &) {
      data.imp_trace.cancel_shift();
      return 0;
    }

    double random_number = rng.preview();
    if (random_number == 0.0) return 0;
    double p_yee = std::abs(1.0 / data.atomic_weight);

    std::tie(new_atomic_weight, new_atomic_reweighting) = data.imp_trace.compute(p_yee, random_number);
    if (new_atomic_weight == 0.0) return 0;

    auto atomic_weight_ratio = new_atomic_weight / data.atomic_weight;
    if (!isfinite(atomic_weight_ratio))
      TRIQS_RUNTIME_ERROR << "(worm shift F) atomic_weight_ratio not finite " << new_atomic_weight << " " << data.atomic_weight;

    return atomic_weight_ratio;
  }

  mc_weight_t move_worm_shift_F::accept() {
    data.imp_trace.confirm_shift();
    if (shift_Q)
      data.worm.tau_Q = tau_new;
    else
      data.worm.tau_cdag = tau_new;
    data.atomic_weight      = new_atomic_weight;
    data.atomic_reweighting = new_atomic_reweighting;
    return 1;
  }

  void move_worm_shift_F::reject() { data.imp_trace.cancel_shift(); }

} // namespace triqs_cthyb
