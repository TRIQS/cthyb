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

#include "./g2_tau.hpp"

namespace cthyb {

  using namespace triqs::gfs;

  measure_g2_tau::measure_g2_tau(int A, int B, g2_tau_view_type g2_tau, qmc_data const &data)
     : data(data), g2_tau(g2_tau), A(A), B(B), average_sign(0.0) {
    g2_tau() = 0.0;
  }

  void measure_g2_tau::accumulate(mc_weight_t sign) {

    sign *= data.atomic_reweighting;
    average_sign += sign;

    foreach (data.dets[A], [&](auto const &j, auto const &i, auto const M_ij) {
      foreach (data.dets[B], [&](auto const &l, auto const &k, auto const M_kl) {

        // lambda for computing a single product term of M_ij and M_kl
        auto compute_M2_product = [&](auto const &i, auto const &j, auto const &k, auto const &l, double sign) {

          double t1 = double(i.first - l.first);
          double t2 = double(j.first - l.first);
          double t3 = double(k.first - l.first);

          // implicit beta-periodicity, but fix the sign properly
          int sign_flips    = int(i.first < l.first) + int(j.first < l.first) + int(k.first < l.first);
          double pre_factor = (sign_flips % 2 ? -sign : sign);

          this->g2_tau[closest_mesh_pt(t1, t2, t3)](i.second, j.second, k.second, l.second) += pre_factor * M_ij * M_kl;
        };

        compute_M2_product(i, j, k, l, sign);
        if (A == B) compute_M2_product(i, l, k, j, -sign);

      })
        ;
    })
      ;
  }

  void measure_g2_tau::collect_results(triqs::mpi::communicator const &comm) {

    average_sign = mpi_all_reduce(average_sign, comm);
    g2_tau       = mpi_all_reduce(g2_tau, comm);

    // Bin volime in imaginary time space
    double dtau0    = std::get<0>(g2_tau.mesh()).delta();
    double dtau1    = std::get<1>(g2_tau.mesh()).delta();
    double dtau2    = std::get<2>(g2_tau.mesh()).delta();
    double dtau_vol = dtau0 * dtau1 * dtau2;

    // Rescale sampled Green's function
    g2_tau = g2_tau / (real(average_sign) * data.config.beta() * dtau_vol);

    // Account for
    // the 1/2 smaller volume of the side bins,
    // the 1/4 smaller volume of the edge bins, and
    // the 1/8 smaller volume of the corner bins.

    auto _ = var_t{};
    int n  = std::get<0>(g2_tau.mesh().components()).size() - 1;

    g2_tau[0][_][_] *= 2.0;
    g2_tau[_][0][_] *= 2.0;
    g2_tau[_][_][0] *= 2.0;

    g2_tau[n][_][_] *= 2.0;
    g2_tau[_][n][_] *= 2.0;
    g2_tau[_][_][n] *= 2.0;
  }
}
