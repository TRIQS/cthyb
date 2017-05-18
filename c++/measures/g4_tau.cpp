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

#include "./g4_tau.hpp"

namespace cthyb {

  using namespace triqs::gfs;

  measure_g4_tau::measure_g4_tau(std::optional<g4_tau_t> &g4_tau_opt, qmc_data const &data, g4_measures_t const &g4_measures)
     : data(data), g4_measures(g4_measures), average_sign(0) {

    double beta = data.config.beta();
    int n_tau = g4_measures.params.measure_g4_n_tau;
    
    gf_mesh<imtime> fermi_tau_mesh{beta, Fermion, n_tau};
    gf_mesh<cartesian_product<imtime, imtime, imtime>> g4_tau_mesh{fermi_tau_mesh, fermi_tau_mesh, fermi_tau_mesh};

    g4_tau_opt = make_block2_gf(g4_tau_mesh, g4_measures.gf_struct);

    g4_tau.rebind(*g4_tau_opt);
    g4_tau() = 0.0;
  }

  void measure_g4_tau::accumulate(mc_weight_t sign) {

    sign *= data.atomic_reweighting;
    average_sign += sign;

    // loop only over block-combinations that should be measured
    for (auto &m : g4_measures()) {

      auto g4_tau_block = g4_tau(m.b1.idx, m.b2.idx);
      
      foreach (data.dets[m.b1.idx], [&](auto const &j, auto const &i, auto const M_ij) {
        foreach (data.dets[m.b2.idx], [&](auto const &l, auto const &k, auto const M_kl) {

          // lambda for computing a single product term of M_ij and M_kl
          auto compute_M2_product = [&](auto const &i, auto const &j, auto const &k, auto const &l, double sign) {

            double t1 = double(i.first - l.first);
            double t2 = double(j.first - l.first);
            double t3 = double(k.first - l.first);

            // implicit beta-periodicity, but fix the sign properly
            int sign_flips    = int(i.first < l.first) + int(j.first < l.first) + int(k.first < l.first);
            double pre_factor = (sign_flips % 2 ? -sign : sign);

            g4_tau_block[closest_mesh_pt(t1, t2, t3)](i.second, j.second, k.second, l.second) += pre_factor * M_ij * M_kl;
          };

          compute_M2_product(i, j, k, l, sign);
          if (m.b1.idx == m.b2.idx) compute_M2_product(i, l, k, j, -sign);

        })
          ;
      })
        ;
    }
  }

  void measure_g4_tau::collect_results(triqs::mpi::communicator const &comm) {

    average_sign = mpi_all_reduce(average_sign, comm);
    g4_tau       = mpi_all_reduce(g4_tau, comm);

    for (auto &g4_tau_block : g4_tau) {
      // Bin volume in imaginary time space
      double dtau0    = std::get<0>(g4_tau_block.mesh()).delta();
      double dtau1    = std::get<1>(g4_tau_block.mesh()).delta();
      double dtau2    = std::get<2>(g4_tau_block.mesh()).delta();
      double dtau_vol = dtau0 * dtau1 * dtau2;

      // Rescale sampled Green's function
      double beta = data.config.beta();
      g4_tau_block = g4_tau_block / (real(average_sign) * beta * dtau_vol);

      // Account for
      // the 1/2 smaller volume of the side bins,
      // the 1/4 smaller volume of the edge bins, and
      // the 1/8 smaller volume of the corner bins.

      auto _ = var_t{};
      int n  = std::get<0>(g4_tau_block.mesh().components()).size() - 1;

      g4_tau_block[0][_][_] *= 2.0;
      g4_tau_block[_][0][_] *= 2.0;
      g4_tau_block[_][_][0] *= 2.0;

      g4_tau_block[n][_][_] *= 2.0;
      g4_tau_block[_][n][_] *= 2.0;
      g4_tau_block[_][_][n] *= 2.0;
    }
  }
}
