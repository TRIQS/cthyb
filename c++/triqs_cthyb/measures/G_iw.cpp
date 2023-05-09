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

#include "./G_iw.hpp"
#include <triqs/mesh/imfreq.hpp>

namespace triqs_cthyb {

  using namespace triqs::gfs;
  using namespace triqs::mesh;

  measure_G_iw ::measure_G_iw(qmc_data const &data, int n_iw, gf_struct_t const &gf_struct, container_set_t &results) : data(data), average_sign(0) {
    results.G_iw_direct = block_gf<imfreq, G_target_t>{{data.config.beta(), Fermion, n_iw}, gf_struct};
    G_iw.rebind(*results.G_iw_direct);
    G_iw() = 0.0;
  }

  void measure_G_iw::accumulate(mc_weight_t s) {
    s *= data.atomic_reweighting;
    average_sign += s;

    for (auto block_idx : range(G_iw.size())) {
      foreach (data.dets[block_idx], [this, s, block_idx](op_t const &x, op_t const &y, det_scalar_t M) {
        // beta-periodicity is implicit in the argument, just fix the sign properly
        auto val      = (y.first >= x.first ? s : -s) * M;
        double dtau   = double(y.first - x.first);
        auto const &m = G_iw[0].mesh();
        long n_iw     = m.last_index();
        auto ex_dt    = std::exp(M_PI / m.domain().beta * dtau * 1i);
        auto ex_2dt   = ex_dt * ex_dt;
        auto ex       = ex_dt;
        for (long n = 0; n < n_iw; ++n, ex *= ex_2dt) { G_iw[block_idx][n](y.second, x.second) += val * ex; }
        //this->G_tau[block_idx][closest_mesh_pt(dtau)](y.second, x.second) += val;
      })
        ;
    }
  }

  void measure_G_iw::collect_results(mpi::communicator const &c) {

    G_iw         = mpi::all_reduce(G_iw, c);
    average_sign = mpi::all_reduce(average_sign, c);
    G_iw /= -real(average_sign);
    // for (auto &gbl : G_iw) {
    //   //double beta = gbl.mesh().domain().beta;
    //   gbl /= -real(average_sign); //  * beta * gbl.mesh().delta();
    // }
  }

} // namespace triqs_cthyb
