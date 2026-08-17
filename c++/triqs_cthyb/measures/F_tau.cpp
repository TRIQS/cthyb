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

#include "./F_tau.hpp"

namespace triqs_cthyb {

  using namespace triqs::gfs;
  using namespace triqs::mesh;

  measure_F_tau::measure_F_tau(qmc_data const &data, int n_tau, gf_struct_t const &gf_struct, container_set_t &results, double worm_eta)
     : data(data), worm_eta(worm_eta), sign_Z(0) {
    results.F_tau_accum = block_gf<imtime, G_target_t>({data.config.beta(), Fermion, n_tau}, gf_struct);
    F_tau.rebind(*results.F_tau_accum);
    F_tau() = 0.0;
  }

  void measure_F_tau::accumulate(mc_weight_t s) {
    s *= data.atomic_reweighting;

    if (data.worm.in_Z()) {
      sign_Z += s;
      return;
    }

    auto const &worm = data.worm;
    double dtau      = double(worm.tau_Q - worm.tau_cdag);
    F_tau[worm.block][closest_mesh_pt(dtau)](worm.inner_Q, worm.inner_cdag) += s;
  }

  void measure_F_tau::collect_results(mpi::communicator const &c) {
    F_tau = mpi::all_reduce(F_tau, c);
    sign_Z = mpi::all_reduce(sign_Z, c);

    for (auto &F_tau_block : F_tau) {
      double beta = F_tau_block.mesh().beta();
      F_tau_block /= -real(sign_Z) * beta * F_tau_block.mesh().delta() * worm_eta;

      int last = F_tau_block.mesh().size() - 1;
      F_tau_block[0] *= 2;
      F_tau_block[last] *= 2;
    }
  }

} // namespace triqs_cthyb
