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

#include "./F_tau_partition.hpp"

#include <map>

namespace triqs_cthyb {

  using namespace triqs::gfs;
  using namespace triqs::mesh;

  measure_F_tau_partition::measure_F_tau_partition(qmc_data const &data, int n_tau, gf_struct_t const &gf_struct, container_set_t &results)
     : data(data), average_sign(0) {
    results.F_tau_partition_accum = block_gf<imtime, G_target_t>({data.config.beta(), Fermion, n_tau}, gf_struct);
    F_tau_partition.rebind(*results.F_tau_partition_accum);
    F_tau_partition() = 0.0;
  }

  void measure_F_tau_partition::accumulate(mc_weight_t s) {
    if (!data.worm.in_Z()) return;

    s *= data.atomic_reweighting;
    average_sign += s;

    auto [w0, rw0] = data.imp_trace.compute();
    auto trace0    = w0 * rw0;

    for (auto block_idx : range(F_tau_partition.size())) {
      auto const &det = data.dets[block_idx];
      long n          = det.size();
      if (n == 0) continue;

      std::map<time_pt, mc_weight_t> trace_ratio;
      for (long j = 0; j < n; ++j) {
        auto y             = det.get_y(j);
        auto const &Q_desc = data.worm.Q_ops[block_idx][y.second];
        if (!Q_desc) continue;

        configuration::oplist_t updated_ops;
        updated_ops.emplace(y.first, *Q_desc);
        try {
          data.imp_trace.try_replace(updated_ops);
          auto [wQ, rwQ]      = data.imp_trace.compute();
          trace_ratio[y.first] = (wQ * rwQ) / trace0;
        } catch (...) {
          data.imp_trace.cancel_replace();
          throw;
        }
        data.imp_trace.cancel_replace();
      }

      foreach (det, [this, s, block_idx, &trace_ratio](op_t const &x, op_t const &y, det_scalar_t M) {
        auto ratio_it = trace_ratio.find(y.first);
        if (ratio_it == trace_ratio.end()) return;

        auto val    = (y.first >= x.first ? s : -s) * M * ratio_it->second;
        double dtau = double(y.first - x.first);
        this->F_tau_partition[block_idx][closest_mesh_pt(dtau)](y.second, x.second) += val;
      });
    }
  }

  void measure_F_tau_partition::collect_results(mpi::communicator const &c) {
    F_tau_partition = mpi::all_reduce(F_tau_partition, c);
    average_sign    = mpi::all_reduce(average_sign, c);

    for (auto &F_tau_block : F_tau_partition) {
      double beta = F_tau_block.mesh().beta();
      F_tau_block /= real(average_sign) * beta * F_tau_block.mesh().delta();

      int last = F_tau_block.mesh().size() - 1;
      F_tau_block[0] *= 2;
      F_tau_block[last] *= 2;
    }
  }

} // namespace triqs_cthyb
