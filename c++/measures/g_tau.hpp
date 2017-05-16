/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, H. U.R. Strand, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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
#include <triqs/gfs.hpp>
#include "./qmc_data.hpp"

namespace cthyb {

  using namespace triqs::gfs;

  // Measure imaginary time Green's function (all blocks)
  struct measure_g_tau {

    g_tau_g_target_t::view_type g_tau;
    qmc_data const &data;
    mc_weight_t average_sign;

    measure_g_tau(std::optional<g_tau_g_target_t> &g_tau_opt, qmc_data const &data, int n_tau, gf_struct_t gf_struct) : data(data), average_sign(0) {
      g_tau_opt = make_block_gf<g_target_t>(gf_mesh<imtime>{data.config.beta(), Fermion, n_tau}, gf_struct);
      g_tau.rebind(*g_tau_opt);
      g_tau() = 0.0;
    }
    // --------------------

    void accumulate(mc_weight_t s) {
      s *= data.atomic_reweighting;
      average_sign += s;

      for (auto block_idx : range(g_tau.size())) {
        foreach (data.dets[block_idx], [this, s, block_idx](op_t const &x, op_t const &y, det_scalar_t M) {
          // beta-periodicity is implicit in the argument, just fix the sign properly
          auto val    = (y.first >= x.first ? s : -s) * M;
          double dtau = double(y.first - x.first);
          this->g_tau[block_idx][closest_mesh_pt(dtau)](y.second, x.second) += val;
        })
          ;
      }
    }
    // ---------------------------------------------

    void collect_results(triqs::mpi::communicator const &c) {

      g_tau        = mpi_all_reduce(g_tau, c);
      average_sign = mpi_all_reduce(average_sign, c);

      for (auto &g_block : g_tau) {
        double beta = g_block.mesh().domain().beta;
        g_block /= -real(average_sign) * beta * g_block.mesh().delta();

        // Multiply first and last bins by 2 to account for full bins
        int last = g_block.mesh().size() - 1;
        g_block[0] *= 2;
        g_block[last] *= 2;

        // Set 1/iw behaviour of tails in G_tau to avoid problems when taking FTs later
        g_block.singularity()(1) = 1.0;
      }
    }
  };
  // ---------------------------------------------
}
