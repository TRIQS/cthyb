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
#pragma once

#include "../qmc_data.hpp"

#include <triqs/utility/exceptions.hpp>

namespace triqs_cthyb::detail {

  struct partition_measurement_context {
    mc_weight_t mc_sign;
    h_scalar_t atomic_weight;
    mc_weight_t average_sign_contribution;
  };

  inline partition_measurement_context make_partition_measurement_context(qmc_data const &data, mc_weight_t mc_sign) {
    auto [atomic_weight, atomic_reweighting] = data.imp_trace.compute();
    if (!isfinite(atomic_weight) || std::abs(atomic_weight) == 0)
      TRIQS_RUNTIME_ERROR << "Partition estimator encountered a zero or non-finite Monte Carlo atomic weight in an accepted Z-sector "
                             "configuration: "
                          << atomic_weight;

    return {mc_sign, atomic_weight, mc_sign * atomic_reweighting};
  }

  inline h_scalar_t replacement_trace_over_atomic_weight(qmc_data const &data, h_scalar_t atomic_weight) {
    // The direct trace-ratio form is numerically unsafe with norm reweighting:
    //
    //   mc_sign * (trace_0 / weight_0) * (trace_Q / trace_0).
    //
    // A sampled configuration can have weight_0 = ||rho_0||_F > 0 but trace_0 = 0,
    // producing 0 * inf or 0 / 0. Cancel trace_0 algebraically first and evaluate
    // the equivalent finite expression mc_sign * trace_Q / weight_0 instead.
    auto [replacement_weight, replacement_reweighting] = data.imp_trace.compute();
    auto replacement_trace                             = replacement_weight * replacement_reweighting;
    auto result                                        = replacement_trace / atomic_weight;

    if (!isfinite(result))
      TRIQS_RUNTIME_ERROR << "Partition estimator replacement trace divided by the Monte Carlo atomic weight is not finite: trace = "
                          << replacement_trace << ", weight = " << atomic_weight;
    return result;
  }

} // namespace triqs_cthyb::detail
