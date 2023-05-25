/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, H. U.R. Strand, M. Ferrero and O. Parcollet
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

#include "./G2_tau.hpp"

namespace triqs_cthyb {

  using namespace triqs::gfs;
  using namespace triqs::mesh;

  measure_G2_tau::measure_G2_tau(std::optional<G2_tau_t> &G2_tau_opt, qmc_data const &data, G2_measures_t const &G2_measures)
     : data(data), average_sign(0), G2_measures(G2_measures) {

    double beta = data.config.beta();

    order     = G2_measures.params.measure_G2_block_order;
    int n_tau = G2_measures.params.measure_G2_n_tau;

    mesh::imtime fermi_tau_mesh{beta, Fermion, n_tau};
    mesh::prod<imtime, imtime, imtime> G2_tau_mesh{fermi_tau_mesh, fermi_tau_mesh, fermi_tau_mesh};

    G2_tau_opt = make_block2_gf(G2_tau_mesh, G2_measures.gf_struct, order);

    G2_tau.rebind(*G2_tau_opt);
    G2_tau() = 0.0;
  }

  void measure_G2_tau::accumulate(mc_weight_t sign) {

    sign *= data.atomic_reweighting;
    average_sign += sign;

    // loop only over block-combinations that should be measured
    for (auto &m : G2_measures()) {

      auto G2_tau_block = G2_tau(m.b1.index(), m.b2.idx);
      bool diag_block   = (m.b1.index() == m.b2.idx);

      foreach (data.dets[m.b1.index()], [&](auto const &i, auto const &j, auto const M_ij) {
        foreach (data.dets[m.b2.index()], [&](auto const &k, auto const &l, auto const M_kl) {

          // lambda for computing a single product term of M_ij and M_kl
          auto compute_M2_product = [&](auto const &i, auto const &j, auto const &k, auto const &l, mc_weight_t sign) {

            double t1 = double(i.first - l.first);
            double t2 = double(j.first - l.first);
            double t3 = double(k.first - l.first);

            // implicit beta-periodicity, but fix the sign properly
            int sign_flips    = int(i.first < l.first) + int(j.first < l.first) + int(k.first < l.first);
            mc_weight_t pre_factor = (sign_flips % 2 ? -sign : sign);

            G2_tau_block[closest_mesh_pt(t1, t2, t3)](i.second, j.second, k.second, l.second) += pre_factor * M_ij * M_kl;
          };

          if (order == block_order::AABB || diag_block) compute_M2_product(i, j, k, l, +sign);
          if (order == block_order::ABBA || diag_block) compute_M2_product(i, l, k, j, -sign);

        })
          ;
      })
        ;
    }
  }

  void measure_G2_tau::collect_results(mpi::communicator const &comm) {

    average_sign = mpi::all_reduce(average_sign, comm);
    G2_tau       = mpi::all_reduce(G2_tau, comm);

    // Rescale sampled Green's function
    double beta = data.config.beta();
    double dtau = std::get<0>(G2_tau(0,0).mesh()).delta();
    G2_tau = G2_tau / (real(average_sign) * beta * std::pow(dtau, 3));

    // Account for
    // the 1/2 smaller volume of the side bins,
    // the 1/4 smaller volume of the edge bins, and
    // the 1/8 smaller volume of the corner bins.
    
    for (auto &G2_tau_block : G2_tau) {
      auto _ = all_t{};
      int n  = std::get<0>(G2_tau_block.mesh().components()).size() - 1;

      G2_tau_block[0, _, _] *= 2.0;
      G2_tau_block[_, 0, _] *= 2.0;
      G2_tau_block[_, _, 0] *= 2.0;

      G2_tau_block[n, _, _] *= 2.0;
      G2_tau_block[_, n, _] *= 2.0;
      G2_tau_block[_, _, n] *= 2.0;
    }
  }

} // namespace triqs_cthyb
