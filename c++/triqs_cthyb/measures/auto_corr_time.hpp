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

  /// Measurement auto-correlation time based on the partition function
  struct measure_auto_corr_time {

    measure_auto_corr_time(qmc_data const &_data, double &_auto_corr_time) : data(_data), auto_corr_time(_auto_corr_time) {}

    void accumulate(mc_weight_t sign) {
      log_accs[0] << sign;
      log_accs[1] << data.config.size();
    }

    void collect_results(mpi::communicator const &comm) {

      constexpr int min_samples = 32;
      auto_corr_time            = 0.0;

      // mean_errors_and_taus all-reduces internally, so every rank obtains the same result.
      for (auto &log_acc : log_accs) {
        // The integrated autocorrelation time is the saturated (largest-bin) tau estimate.
        auto [mean, errs, taus, effs] = log_acc.mean_errors_and_taus(comm, min_samples);
        // tau is NaN when an observable has zero variance (e.g. sign in sign-problem-free runs).
        if (!taus.empty() && std::isfinite(taus.back())) auto_corr_time = std::max(auto_corr_time, taus.back());
      }
    }

    private:
    qmc_data const &data;
    double &auto_corr_time;

    // Initialize one complex log accumulator for each observable to use for the autocorrelation analysis
    std::vector<triqs::stat::log_binning<dcomplex>> log_accs = {2, {dcomplex{0.0}, -1}};
  };

} // namespace triqs_cthyb
