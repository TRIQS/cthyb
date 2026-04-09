/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2021-2025, Simons Foundation
 *    authors: N. Wentzell
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

#include "./densities.hpp"

#include <cmath>

namespace triqs_cthyb {

  measure_densities::measure_densities(qmc_data &data, gf_struct_t gf_struct, bool measure_densities_flag,
                                       double &auto_corr_time, bool &auto_corr_time_converged,
                                       std::optional<std::map<std::string, nda::array<double, 1>>> &densities,
                                       std::optional<std::map<std::string, nda::array<double, 1>>> &densities_errors)
     : data(data),
       gf_struct(std::move(gf_struct)),
       measure_densities_(measure_densities_flag),
       auto_corr_time(auto_corr_time),
       auto_corr_time_converged(auto_corr_time_converged),
       densities_(densities),
       densities_errors_(densities_errors) {

    // Log-binning for auto-correlation: [0] = perturbation order, [1..n_blocks] = per-block det size
    int n_blocks = this->gf_struct.size();
    log_accs_.reserve(1 + n_blocks);
    for (int i = 0; i < 1 + n_blocks; ++i) log_accs_.emplace_back(dcomplex{0.0}, -1);

    // Attach n_a = c†_a c_a auxiliary operators and store their indices
    if (measure_densities_) {
      using triqs::operators::n;
      for (auto const &[bl_name, bl_size] : this->gf_struct) {
        dens_bins_.emplace_back(nda::zeros<dcomplex>(bl_size), 128, 1);
        for (int a = 0; a < bl_size; ++a) {
          auto op_d = data.imp_trace.attach_aux_operator(n(bl_name, a));
          n_op_indices_.push_back(-op_d.linear_index - 1);
        }
      }
    }
  }

  void measure_densities::accumulate(mc_weight_t sign) {

    // Log-bin the total perturbation order
    log_accs_[0] << double(data.config.size() / 2);

    // Log-bin per-block det sizes
    for (int b = 0; b < static_cast<int>(data.dets.size()); ++b) { log_accs_[1 + b] << double(data.dets[b].size()); }

    if (!measure_densities_) return;

    Z += sign;
    ++N_;

    // Compute bare trace (caches root matrices for trace_with_aux_op)
    auto [bare_w, bare_rw] = data.imp_trace.compute();
    auto Z_cfg              = bare_w * bare_rw;

    // Measure per-orbital densities from cached root matrices — no tree modifications
    int op_idx = 0;
    for (int b = 0; b < static_cast<int>(gf_struct.size()); ++b) {
      int bl_size = gf_struct[b].second;
      auto step   = nda::array<dcomplex, 1>(bl_size);

      for (int a = 0; a < bl_size; ++a) {
        step(a) = sign * data.imp_trace.trace_with_aux_op(n_op_indices_[op_idx]) / Z_cfg;
        ++op_idx;
      }
      dens_bins_[b] << step;
    }
  }

  void measure_densities::collect_results(mpi::communicator const &comm) {
    using triqs::stat::log_binning;

    // Auto-correlation time from log-binning (always active).
    // tau saturates at large bin size. We report it at the deepest bin level with >= min_samples effective
    // samples and flag it as a lower bound if its finer neighbour (half the bin size) is more than n_sigma
    // errors below it. The error of tau from M samples is d_tau = (tau + 1/2) * sqrt(2 / (M - 1)) (variance-
    // of-a-variance); the deeper level dominates the uncertainty, so we use its d_tau as the yardstick.
    constexpr int min_samples = 64;  // report level: ~9% error on the error bar
    constexpr double n_sigma  = 2.0; // flag as "rising" if tau grows by > n_sigma errors per bin doubling (~95%)

    auto d_tau = [](double tau, long M) { return (tau + 0.5) * std::sqrt(2.0 / static_cast<double>(M - 1)); };

    // mean_errors_and_taus all-reduces internally, so every rank obtains the same result.
    auto_corr_time           = 0.0;
    auto_corr_time_converged = true;

    for (auto &log_acc : log_accs_) {
      auto [mean, errs, taus, effs] = log_acc.mean_errors_and_taus(comm, min_samples);
      log_acc                       = log_binning<dcomplex>{dcomplex{0.0}, -1};

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
      if (tau_a - taus[n - 2] > n_sigma * d_tau(tau_a, effs[n - 1])) auto_corr_time_converged = false; // still rising -> lower bound
    }

    if (!measure_densities_) return;

    Z  = mpi::all_reduce(Z, comm);
    N_ = mpi::all_reduce(N_, comm);

    // Compute densities and error bars from linear binning
    std::map<std::string, nda::array<double, 1>> densities;
    std::map<std::string, nda::array<double, 1>> densities_errors;
    auto norm = std::abs(Z / double(N_));

    for (int b = 0; auto const &[bl_name, bl_size] : gf_struct) {
      auto [m, err, tau] = dens_bins_[b].mean_error_and_tau(comm);

      // Mean density = mean * N / Z (convert from sign-weighted per-step to normalized)
      densities[bl_name] = nda::array<double, 1>(nda::real(m / norm));
      densities_errors[bl_name] = nda::array<double, 1>(nda::abs(err) / norm);
      ++b;
    }

    if (comm.rank() == 0) {
      std::cout << "Densities:" << std::endl;
      for (auto &[bl, dens] : densities) std::cout << "  " << bl << ": " << dens << std::endl;
      std::cout << "Densities errors:" << std::endl;
      for (auto &[bl, err] : densities_errors) std::cout << "  " << bl << ": " << err << std::endl;
    }

    densities_        = std::move(densities);
    densities_errors_ = std::move(densities_errors);
  }

} // namespace triqs_cthyb
