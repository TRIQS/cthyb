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
#include <triqs/utility/legendre.hpp>
#include "qmc_data.hpp"
#include <boost/serialization/complex.hpp>
#include <boost/mpi/collectives.hpp>

namespace cthyb {

using namespace triqs::gfs;

// Measure Legendre Green's function (one block)
struct measure_g_legendre {
 using mc_sign_type = std::complex<double>;

 qmc_data const& data;
 gf_view<legendre> g_l;
 int a_level;
 double beta;
 mc_sign_type z;
 int64_t num;
 mc_sign_type average_sign;

 measure_g_legendre(int a_level, gf_view<legendre> g_l, qmc_data const& data)
    : data(data), g_l(g_l), a_level(a_level), beta(data.config.beta()) {
  z = 0;
  num = 0;
 }
 // --------------------

 void accumulate(mc_sign_type s) {
  num += 1;
  if (num < 0) TRIQS_RUNTIME_ERROR << " Overflow of counter ";

  auto corr = real(this->data.atomic_corr.full_trace_over_estimator());
  if (!std::isfinite(corr)) TRIQS_RUNTIME_ERROR << " measure g :corr not finite" << corr;

  z += s * corr;

  auto Tn = triqs::utility::legendre_generator();
  
  foreach(data.dets[a_level], [this, corr, s, &Tn](std::pair<time_pt, int> const& x, std::pair<time_pt, int> const& y, double M) {
   double poly_arg = 2*double(y.first - x.first)/beta - 1.0;   
   Tn.reset(poly_arg);
   
   double val = (y.first >= x.first ? real(s) : -real(s)) * M * corr;
   for (auto l : g_l.mesh()) this->g_l[l](y.second, x.second) += val*Tn.next();
  });
 }
 // ---------------------------------------------

 void collect_results(boost::mpi::communicator const& c) {

  boost::mpi::all_reduce(c, z, z, std::c14::plus<>());
  boost::mpi::all_reduce(c, num, num, std::c14::plus<>());
  average_sign = z / num;
    
  auto g_l_copy = make_clone(g_l);
  boost::mpi::all_reduce(c, g_l_copy, g_l_copy, std::c14::plus<>());
  for (auto l : g_l_copy.mesh()) g_l[l] = -(sqrt(2.0*l+1.0)/(real(z)*beta)) * g_l_copy[l];
 }
};
}
