/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2021, Simons Foundation
 *    author: N. Wentzell
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

#include <triqs/stat/log_binning.hpp>

#include "../qmc_data.hpp"

#include <cmath>

namespace triqs_cthyb {

  /// Auto-correlation time from the partition function (one sample per cycle, so tau is in units of MC cycles)
  struct measure_auto_corr_time {

    measure_auto_corr_time(qmc_data const &_data, double &_auto_corr_time, bool &_auto_corr_time_converged)
       : data(_data), auto_corr_time(_auto_corr_time), auto_corr_time_converged(_auto_corr_time_converged) {}

    void accumulate(mc_weight_t sign) {
      log_accs[0] << sign;
      log_accs[1] << data.config.size();
    }

    void collect_results(mpi::communicator const &comm) {

      // tau saturates at large bin size. We report it at the deepest bin level with >= min_samples effective
      // samples and flag it as a lower bound if its finer neighbour (half the bin size) is more than n_sigma
      // errors below it. The error of tau from M samples is d_tau = (tau + 1/2) * sqrt(2 / (M - 1)) (variance-
      // of-a-variance); the deeper level dominates the uncertainty, so we use its d_tau as the yardstick.
      constexpr int min_samples = 64;  // report level: ~9% error on the error bar
      constexpr double n_sigma  = 2.0; // flag as "rising" if tau grows by > n_sigma errors per bin doubling (~95%)

      // mean_errors_and_taus all-reduces internally, so every rank obtains the same result.
      auto_corr_time           = 0.0;
      auto_corr_time_converged = true;

      auto d_tau = [](double tau, long M) { return (tau + 0.5) * std::sqrt(2.0 / static_cast<double>(M - 1)); };

      for (auto &log_acc : log_accs) {
        auto [mean, errs, taus, effs] = log_acc.mean_errors_and_taus(comm, min_samples);

        // No bin level reached min_samples effective samples: the run is too short to estimate tau at all,
        // so we cannot claim it has saturated.
        if (taus.empty()) {
          auto_corr_time_converged = false;
          continue;
        }

        // Zero variance (e.g. the sign in a sign-problem-free run): tau is genuinely ~0, treat as saturated.
        if (!std::isfinite(taus.back())) continue;

        double const tau_a = taus.back();
        auto_corr_time     = std::max(auto_corr_time, tau_a);

        // Compare against the finer neighbour to check whether tau has saturated.
        auto const n = taus.size();
        if (n < 2 || !std::isfinite(taus[n - 2])) {
          auto_corr_time_converged = false; // no neighbour to confirm saturation
          continue;
        }
        double const tau_b = taus[n - 2];
        if (tau_a - tau_b > n_sigma * d_tau(tau_a, effs[n - 1])) auto_corr_time_converged = false; // still rising -> lower bound
      }
    }

    private:
    qmc_data const &data;
    double &auto_corr_time;
    bool &auto_corr_time_converged;

    // One complex log-binning accumulator per observable (partition-function sign and perturbation order).
    std::vector<triqs::stat::log_binning<dcomplex>> log_accs = {2, {dcomplex{0.0}, -1}};
  };

} // namespace triqs_cthyb
