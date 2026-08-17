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

#include "./F_l_partition.hpp"

#include <map>

namespace triqs_cthyb {

  using namespace triqs::gfs;
  using namespace triqs::mesh;

  measure_F_l_partition::measure_F_l_partition(std::optional<G_l_t> &F_l_partition_opt, qmc_data const &data, int n_l,
                                               gf_struct_t const &gf_struct)
     : data(data), average_sign(0) {
    F_l_partition_opt = block_gf<legendre>{{data.config.beta(), Fermion, n_l}, gf_struct};
    F_l_partition.rebind(*F_l_partition_opt);
    F_l_partition() = 0.0;
  }

  void measure_F_l_partition::accumulate(mc_weight_t s) {
    if (!data.worm.in_Z()) return;

    s *= data.atomic_reweighting;
    average_sign += s;

    auto [w0, rw0] = data.imp_trace.compute();
    auto trace0    = w0 * rw0;
    double beta    = data.config.beta();

    for (auto block_idx : range(F_l_partition.size())) {
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

      foreach (det, [this, s, block_idx, beta, &trace_ratio](op_t const &x, op_t const &y, det_scalar_t M) {
        auto ratio_it = trace_ratio.find(y.first);
        if (ratio_it == trace_ratio.end()) return;

        double poly_arg = 2 * double(y.first - x.first) / beta - 1.0;
        auto Tn         = triqs::utility::legendre_generator();
        Tn.reset(poly_arg);

        auto val = (y.first >= x.first ? s : -s) * M * ratio_it->second;
        for (auto l : this->F_l_partition[block_idx].mesh())
          this->F_l_partition[block_idx][l](y.second, x.second) += val * Tn.next();
      });
    }
  }

  void measure_F_l_partition::collect_results(mpi::communicator const &c) {
    F_l_partition = mpi::all_reduce(F_l_partition, c);
    average_sign  = mpi::all_reduce(average_sign, c);

    double beta = data.config.beta();
    for (auto &F_l_block : F_l_partition)
      for (auto l : F_l_block.mesh()) F_l_block[l] *= sqrt(2.0 * l.index() + 1.0) / (real(average_sign) * beta);
  }

} // namespace triqs_cthyb
