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
#include <cmath>

namespace cthyb {

  using namespace triqs::gfs;

  measure_g2_tau::measure_g2_tau(int A, int B, g2_view_type g2, qmc_data const &data) : data(data), g2(g2), A(A), B(B), diag_block(A == B) {

    g2() = 0.0;
    z    = 0;
    num  = 0;
  }

  void measure_g2_tau::accumulate(mc_weight_t s) {
    num += 1;
    if (num < 0) TRIQS_RUNTIME_ERROR << " Overflow of counter ";

    s *= data.atomic_reweighting;
    z += s;

    auto const &det_A = data.dets[A];
    auto const &det_B = data.dets[B];
    if (det_A.size() == 0 || det_B.size() == 0) return;

    using idx_t = std::pair<time_pt, int>;

    foreach (data.dets[A], [&](idx_t const &j, idx_t const &i, det_scalar_t const M_ij) {
      foreach (data.dets[B], [&](idx_t const &l, idx_t const &k, det_scalar_t const M_kl) {

        auto add = [&](idx_t const &i, idx_t const &j, idx_t const &k, idx_t const &l, double sign) {

          double t1 = double(i.first - l.first);
          double t2 = double(j.first - l.first);
          double t3 = double(k.first - l.first);

          // implicit beta-periodicity, but fix the sign properly
          int sign_flips = int(i.first < l.first) + int(j.first < l.first) + int(k.first < l.first);
          double factor  = (sign_flips % 2 ? -sign : sign);

          this->g2[closest_mesh_pt(t1, t2, t3)](i.second, j.second, k.second, l.second) += factor * M_ij * M_kl;
        };


        add(i, j, k, l, s);
        if (A == B) add(i, l, k, j, -s);

      })
        ;
    })
      ;
  }

  void measure_g2_tau::collect_results(triqs::mpi::communicator const &c) {

    z = mpi_all_reduce(z, c);
    g2 = mpi_all_reduce(g2, c);

    // Fixme: use product reduction on delta()
    double dtau0 = std::get<0>(g2.mesh().components()).delta();
    double dtau1 = std::get<1>(g2.mesh().components()).delta();
    double dtau2 = std::get<2>(g2.mesh().components()).delta();
    double dtau3 = dtau0 * dtau1 * dtau2;

    g2 = g2 / (real(z) * data.config.beta() * dtau3); // Q: tau mesh delta tau?

    // Account for edge bins beeing smaller
    auto _ = var_t{};
    int n  = std::get<0>(g2.mesh().components()).size() - 1;

    // There are 6 sides of the cube
    
    g2[0][_][_] *= 2.0;
    g2[_][0][_] *= 2.0;
    g2[_][_][0] *= 2.0;

    g2[n][_][_] *= 2.0;
    g2[_][n][_] *= 2.0;
    g2[_][_][n] *= 2.0;

    /*
    // There are 3*4 = 12 edges of the cube

    g2[0][0][_] *= 2.0;
    g2[_][0][0] *= 2.0;
    g2[0][_][0] *= 2.0;

    g2[n][n][_] *= 2.0;
    g2[_][n][n] *= 2.0;
    g2[n][_][n] *= 2.0;

    g2[n][0][_] *= 2.0;
    g2[_][n][0] *= 2.0;
    g2[n][_][0] *= 2.0;

    g2[0][n][_] *= 2.0;
    g2[_][0][n] *= 2.0;
    g2[0][_][n] *= 2.0;

    // there are 8 corners of the cube

    g2[0][0][0] *= 2.0; // 1
    g2[n][0][0] *= 2.0; // 2
    g2[0][n][0] *= 2.0; // 3
    g2[0][0][n] *= 2.0; // 4
    g2[n][n][0] *= 2.0; // 5
    g2[n][0][n] *= 2.0; // 6
    g2[0][n][n] *= 2.0; // 7
    g2[n][n][n] *= 2.0; // 8
    */
    
  }  
}
