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

#include <cstdint>

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/utility/legendre.hpp>

#include "../container_set.hpp"
#include "../qmc_data.hpp"

namespace triqs_cthyb {

  using namespace triqs::gfs;
  using namespace triqs::mesh;

  namespace detail {

    /// Per-rank cadence over all accumulation events, including events in non-Z sectors.
    class partition_measurement_schedule {
      public:
      explicit partition_measurement_schedule(long stride) : stride(stride) {
        if (stride < 1) TRIQS_RUNTIME_ERROR << "measure_F_partition_stride must be at least 1, got " << stride;
      }

      bool select_next() { return event_index++ % stride == 0; }

      private:
      std::uint64_t event_index = 0;
      std::uint64_t stride;
    };

  } // namespace detail

  /// Fused Z-sector hybridization-line replacement estimator for optional imaginary-time and Legendre outputs.
  class measure_F_partition {
    public:
    measure_F_partition(qmc_data const &data, int n_tau, int n_l, gf_struct_t const &gf_struct, container_set_t &results, bool measure_tau,
                        bool measure_l, long stride);
    void accumulate(mc_weight_t s);
    void collect_results(mpi::communicator const &c);

    private:
    qmc_data const &data;
    G_tau_G_target_t *F_tau_partition = nullptr;
    G_l_t *F_l_partition              = nullptr;
    mc_weight_t Z_normalization       = 0;
    double absolute_normalization     = 0;
    long selected_Z_events            = 0;
    detail::partition_measurement_schedule schedule;
  };

} // namespace triqs_cthyb
