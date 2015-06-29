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
#include "qmc_data.hpp"
#include <boost/serialization/complex.hpp>
#include <boost/mpi/collectives.hpp>

namespace cthyb {

using namespace triqs::gfs;

// Measure imaginary time Green's function (one block)
struct measure_g {
 using mc_sign_type = std::complex<double>;

 qmc_data const& data;
 gf_view<imtime> g_tau;
 int a_level;
 double beta;
 mc_sign_type z;
 int64_t num;
 mc_sign_type average_sign;

 measure_g(int a_level, gf_view<imtime> g_tau, qmc_data const& data)
    : data(data), g_tau(g_tau), a_level(a_level), beta(data.config.beta()) {
  g_tau() = 0.0;
  z = 0;
  num = 0;
 }
 // --------------------

 void accumulate(mc_sign_type s) {
  num += 1;
  if (num < 0) TRIQS_RUNTIME_ERROR << " Overflow of counter ";

  auto corr = real(this->data.imp_trace.full_trace_over_estimator());
  if (!std::isfinite(corr)) TRIQS_RUNTIME_ERROR << " measure g :corr not finite" << corr;

  z += s * corr;

  foreach(data.dets[a_level], [this, corr, s](std::pair<time_pt, int> const& x, std::pair<time_pt, int> const& y, double M) {
   // beta-periodicity is implicit in the argument, just fix the sign properly
   this->g_tau[closest_mesh_pt(double(y.first - x.first))](y.second, x.second) +=
       (y.first >= x.first ? real(s) : -real(s)) * M * corr;
  });
 }
 // ---------------------------------------------

 void collect_results(boost::mpi::communicator const& c) {

  int64_t total_num;
  mc_sign_type total_z;
  boost::mpi::all_reduce(c, z, total_z, std::c14::plus<>());
  boost::mpi::all_reduce(c, num, total_num, std::c14::plus<>());
  average_sign = total_z / total_num;
  // Multiply first and last bins by 2 to account for full bins
  g_tau[0] = g_tau[0] * 2;
  g_tau[g_tau.mesh().size()-1] = g_tau[g_tau.mesh().size()-1] * 2;
  // Need copies, because all_reduce wants default-constructible types
  auto g_tau_in = make_clone(g_tau);
  auto g_tau_out = make_clone(g_tau);
  boost::mpi::all_reduce(c, g_tau_in, g_tau_out, std::c14::plus<>());
  g_tau = g_tau_out / (-real(total_z) * data.config.beta() * g_tau_out.mesh().delta());
  // Set 1/iw behaviour of tails in G_tau to avoid problems when taking FTs later
  g_tau.singularity()(1) = 1.0;
 }
};
}
