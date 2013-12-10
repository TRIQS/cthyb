/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by I. Krivenko, M. Ferrero, O. Parcollet
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
#include <boost/mpi/collectives.hpp>

namespace cthyb_krylov {

using namespace triqs::gfs;

// Measure imaginary time Green's function (one block)
struct measure_g {

 using mc_sign_type = std::complex<double>;

 qmc_data const& data;
 gf_view<imtime> g_tau;
 int a_level;
 double beta;

 mc_sign_type z;
 long long num;
 mc_sign_type average_sign;

 measure_g(int a_level, gf_view<imtime> g_tau, qmc_data const& data)
    : data(data), g_tau(g_tau), a_level(a_level), beta(data.config.beta()) {
  z = 0;
  num = 0;
 }

 // ---- TO ADJUST with lambda ....
 
 void accumulate(mc_sign_type s) {
  z += s;
  num += 1;

  auto const& tau_mesh = g_tau.mesh();
  auto const& det = data.dets[a_level];

  int size = det.size();
  for (int i = 0; i < size; ++i) {
   auto x_i = det.get_x(i);
   for (int j = 0; j < size; ++j) {
    auto y_j = det.get_y(j);
    auto M = det.inverse_matrix(j, i);
    if (y_j.first >= x_i.first)
     g_tau[tau_mesh.nearest_index(double(y_j.first - x_i.first))](y_j.second, x_i.second) += real(s) * M;
    else
     // beta-periodicity is implicit, just fix the sign
     g_tau[tau_mesh.nearest_index(double(y_j.first - x_i.first))](y_j.second, x_i.second) -= real(s) * M;
   }
  }
 }

 // ---------------------------------------------
 
 void collect_results(boost::mpi::communicator const& c) {

  boost::mpi::all_reduce(c, z, z, std::plus<mc_sign_type>());
  boost::mpi::all_reduce(c, num, num, std::plus<long long>());
  average_sign = z / num;

  // Need a copy, because all_reduce wants default-constructible types
  auto g_tau_copy = triqs::make_clone(g_tau);
  boost::mpi::all_reduce(c, g_tau_copy, g_tau_copy, std::plus<triqs::gfs::gf<imtime>>());
  g_tau = g_tau_copy / (-real(z) * data.config.beta() * g_tau_copy.mesh().delta());
 }
};
}
