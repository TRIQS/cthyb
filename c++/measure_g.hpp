/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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

// Measure imaginary time Green's function (one block)
struct measure_g {

 qmc_data const& data;
 gf_view<imtime, g_target_t> g_tau;
 int a_level;
 mc_weight_t z;
 int64_t num;
 mc_weight_t average_sign;

 measure_g(int a_level, gf_view<imtime, g_target_t> g_tau, qmc_data const& data)
    : data(data), g_tau(g_tau), a_level(a_level) {
  g_tau() = 0.0;
  z = 0;
  num = 0;
 }
 // --------------------

 void accumulate(mc_weight_t s) {
  num += 1;
  if (num < 0) TRIQS_RUNTIME_ERROR << " Overflow of counter ";

  s *= data.atomic_reweighting;
  z += s;

  foreach(data.dets[a_level], [this, s](std::pair<time_pt, int> const& x, std::pair<time_pt, int> const& y, det_scalar_t M) {
   // beta-periodicity is implicit in the argument, just fix the sign properly
   this->g_tau[closest_mesh_pt(double(y.first - x.first))](y.second, x.second) +=
       (y.first >= x.first ? s : -s) * M;
  });
 }
 // ---------------------------------------------

 void collect_results(triqs::mpi::communicator const& c) {

  z = mpi_all_reduce(z,c);
  // Multiply first and last bins by 2 to account for full bins
  g_tau[0] = g_tau[0] * 2;
  g_tau[g_tau.mesh().size() - 1] = g_tau[g_tau.mesh().size() - 1] * 2;
  g_tau = mpi_all_reduce(g_tau, c);
  g_tau = g_tau / (-real(z) * data.config.beta() * g_tau.mesh().delta());
  // Set 1/iw behaviour of tails in G_tau to avoid problems when taking FTs later
  g_tau.singularity()(1) = 1.0;
 }
};
}
