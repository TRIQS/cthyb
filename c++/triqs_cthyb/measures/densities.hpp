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
#pragma once

#include <triqs/stat/log_binning.hpp>
#include <triqs/stat/lin_binning.hpp>

#include "../qmc_data.hpp"

namespace triqs_cthyb {

  /// Unified densities measure: auto-correlation time + optional per-orbital densities with error bars
  struct measure_densities {

    measure_densities(qmc_data &data, gf_struct_t gf_struct, bool measure_densities_flag, double &auto_corr_time,
                      bool &auto_corr_time_converged, std::optional<std::map<std::string, nda::array<double, 1>>> &densities,
                      std::optional<std::map<std::string, nda::array<double, 1>>> &densities_errors);

    void accumulate(mc_weight_t sign);
    void collect_results(mpi::communicator const &comm);

    private:
    qmc_data &data;
    gf_struct_t gf_struct;
    bool measure_densities_;

    double &auto_corr_time;
    bool &auto_corr_time_converged;
    std::optional<std::map<std::string, nda::array<double, 1>>> &densities_;
    std::optional<std::map<std::string, nda::array<double, 1>>> &densities_errors_;

    mc_weight_t Z = 0;
    long N_       = 0;

    // Log-binning for auto-correlation: [0] = perturbation order, [1..n_blocks] = per-block det size
    std::vector<triqs::stat::log_binning<dcomplex>> log_accs_;

    // Auxiliary operator indices for n_a = c†_a c_a (one per orbital, flattened across blocks)
    std::vector<int> n_op_indices_;

    // Linear binning for density errors (one per block)
    std::vector<triqs::stat::lin_binning<nda::array<dcomplex, 1>>> dens_bins_;
  };

} // namespace triqs_cthyb
