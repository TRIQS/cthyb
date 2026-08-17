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

#include "./F_partition.hpp"
#include "./F_partition_common.hpp"

#include <type_traits>

namespace triqs_cthyb {

  using namespace triqs::gfs;
  using namespace triqs::mesh;

  measure_F_partition::measure_F_partition(qmc_data const &data, int n_tau, int n_l, gf_struct_t const &gf_struct, container_set_t &results,
                                           bool measure_tau, bool measure_l, long stride)
     : data(data), schedule(stride) {
    if (!measure_tau && !measure_l) TRIQS_RUNTIME_ERROR << "measure_F_partition requires at least one output";

    if (measure_tau) {
      results.F_tau_partition_accum = block_gf<imtime, G_target_t>({data.config.beta(), Fermion, n_tau}, gf_struct);
      F_tau_partition               = &*results.F_tau_partition_accum;
      (*F_tau_partition)()          = 0.0;
    }
    if (measure_l) {
      results.F_l_partition = block_gf<legendre>{{data.config.beta(), Fermion, n_l}, gf_struct};
      F_l_partition         = &*results.F_l_partition;
      (*F_l_partition)()    = 0.0;
    }
  }

  void measure_F_partition::accumulate(mc_weight_t s) {
    if (!schedule.select_next() || !data.worm.in_Z()) return;

    auto const context = detail::make_partition_measurement_context(data, s);
    Z_normalization += context.average_sign_contribution;
    absolute_normalization += std::abs(data.atomic_reweighting);
    ++selected_Z_events;

    double const beta          = data.config.beta();
    std::vector<impurity_trace::single_replacement_request> requests;
    for (auto block_idx : range(data.dets.size())) {
      auto const &det = data.dets[block_idx];
      auto const &ys  = det.get_y_internal_order();
      for (long j = 0; j < det.size(); ++j) {
        auto const &y      = ys[j];
        auto const &Q_desc = data.worm.Q_ops[block_idx][y.second];
        if (!Q_desc) continue;
        requests.push_back({y.first, *Q_desc});
      }
    }
    auto const replacement_ratios = detail::replacement_traces_over_atomic_weight(data, requests, context.atomic_weight);

    auto accumulate_outputs = [&](auto tau_tag, auto legendre_tag) {
      constexpr bool measure_tau      = decltype(tau_tag)::value;
      constexpr bool measure_legendre = decltype(legendre_tag)::value;
      using legendre_generator_t       = decltype(triqs::utility::legendre_generator());
      std::optional<legendre_generator_t> Tn;
      if constexpr (measure_legendre) Tn.emplace(triqs::utility::legendre_generator());

      std::size_t request_idx = 0;
      for (auto block_idx : range(data.dets.size())) {
        auto const &det = data.dets[block_idx];
        auto const &xs  = det.get_x_internal_order();
        auto const &ys  = det.get_y_internal_order();
        auto const M    = det.inverse_matrix_internal_order();

        for (long j = 0; j < det.size(); ++j) {
          auto const &y      = ys[j];
          auto const &Q_desc = data.worm.Q_ops[block_idx][y.second];
          if (!Q_desc) continue;
          auto const replacement_ratio = replacement_ratios[request_idx++];

          // y-major traversal reuses one replacement trace for all x endpoints and follows
          // det_manip's internal-order convention M(j, i), identical to foreach(det, ...).
          for (long i = 0; i < det.size(); ++i) {
            auto const &x  = xs[i];
            auto const val = (y.first >= x.first ? context.mc_sign : -context.mc_sign) * M(j, i) * replacement_ratio;

            if constexpr (measure_tau) {
              double const dtau = double(y.first - x.first);
              auto &F_tau_block   = (*F_tau_partition)[block_idx];
              auto const tau_idx  = F_tau_block.mesh().to_data_index(closest_mesh_pt(dtau));
              F_tau_block.data()(tau_idx, y.second, x.second) += val;
            }

            if constexpr (measure_legendre) {
              double const poly_arg = 2 * double(y.first - x.first) / beta - 1.0;
              Tn->reset(poly_arg);
              auto &F_l_block = (*F_l_partition)[block_idx];
              for (long l = 0; l < F_l_block.mesh().size(); ++l) F_l_block.data()(l, y.second, x.second) += val * Tn->next();
            }
          }
        }
      }
    };

    if (F_tau_partition && F_l_partition)
      accumulate_outputs(std::true_type{}, std::true_type{});
    else if (F_tau_partition)
      accumulate_outputs(std::true_type{}, std::false_type{});
    else
      accumulate_outputs(std::false_type{}, std::true_type{});
  }

  void measure_F_partition::collect_results(mpi::communicator const &c) {
    if (F_tau_partition) *F_tau_partition = mpi::all_reduce(*F_tau_partition, c);
    if (F_l_partition) *F_l_partition = mpi::all_reduce(*F_l_partition, c);
    Z_normalization        = mpi::all_reduce(Z_normalization, c);
    absolute_normalization = mpi::all_reduce(absolute_normalization, c);
    selected_Z_events      = mpi::all_reduce(selected_Z_events, c);

    auto const denominator = real(Z_normalization);
    if (selected_Z_events == 0 || absolute_normalization == 0) TRIQS_RUNTIME_ERROR << "Partition estimator collected no selected Z-sector samples";
    if (!isfinite(denominator) || denominator == 0)
      TRIQS_RUNTIME_ERROR << "Partition estimator selected-sample normalization is zero or non-finite: " << denominator << " (absolute normalization "
                          << absolute_normalization << ")";

    double const beta = data.config.beta();
    if (F_tau_partition) {
      for (auto &F_tau_block : *F_tau_partition) {
        F_tau_block /= denominator * beta * F_tau_block.mesh().delta();

        int const last = F_tau_block.mesh().size() - 1;
        F_tau_block[0] *= 2;
        F_tau_block[last] *= 2;
      }
    }

    if (F_l_partition)
      for (auto &F_l_block : *F_l_partition)
        for (auto l : F_l_block.mesh()) F_l_block[l] *= sqrt(2.0 * l.index() + 1.0) / (denominator * beta);
  }

} // namespace triqs_cthyb
