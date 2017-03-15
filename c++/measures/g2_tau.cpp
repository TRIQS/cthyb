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

    foreach (data.dets[A], [this, s](std::pair<time_pt, int> const &i, std::pair<time_pt, int> const &j, det_scalar_t M_ij) {
      foreach (data.dets[B], [this, s, &i, &j, &M_ij](std::pair<time_pt, int> const &k, std::pair<time_pt, int> const &l, det_scalar_t M_kl) {

        // can we do without using closest_mesh_pt()?
        auto sign = [](double tau) { return (tau > 0. ? 1. : -1.); };

        // Direct term
        {
          double t1 = double(i.first - l.first);
          double t2 = double(j.first - l.first);
          double t3 = double(k.first - l.first);

          // implicit beta-periodicity, but fix the sign properly
          double tau_sign = sign(t1) * sign(t2) * sign(t3);

          this->g2[closest_mesh_pt(t1, t2, t3)](i.second, j.second, k.second, l.second) += tau_sign * s * M_ij * M_kl;
        }
        // TODO: ADD SECOND TERM IN THE G2 SAMPLING (M_il, M_kj)

        // Exchange term (with extra factor -1)
        {
          double t1 = double(i.first - j.first);
          double t2 = double(l.first - j.first);
          double t3 = double(k.first - j.first);

          // implicit beta-periodicity, but fix the sign properly
          double tau_sign = sign(t1) * sign(t2) * sign(t3);

          this->g2[closest_mesh_pt(t1, t2, t3)](i.second, l.second, k.second, j.second) += -tau_sign * s * M_ij * M_kl;
        }

      })
        ;
    })
      ;
  }

  void measure_g2_tau::collect_results(triqs::mpi::communicator const &c) {

    z = mpi_all_reduce(z, c);

    // Fixme! Account for bins at the boundaries

    g2 = mpi_all_reduce(g2, c);
    // g2 = g2 / (-real(z) * data.config.beta() * g2.mesh().delta());

    // Fixme: use product reduction on delta()
    double dtau0 = std::get<0>(g2.mesh().components()).delta();
    double dtau1 = std::get<1>(g2.mesh().components()).delta();
    double dtau2 = std::get<2>(g2.mesh().components()).delta();
    double dtau3 = dtau0 * dtau1 * dtau2;

    g2 = g2 / (-real(z) * data.config.beta() * dtau3); // Q: tau mesh delta tau?
  }
}
