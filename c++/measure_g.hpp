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
 int64_t num;
 mc_sign_type average_sign;

 measure_g(int a_level, gf_view<imtime> g_tau, qmc_data const& data)
    : data(data), g_tau(g_tau), a_level(a_level), beta(data.config.beta()) {
  z = 0;
  num = 0;
 }
 // --------------------

 void accumulate(mc_sign_type s) {
  z += s;
  num += 1;
  if (num < 0) TRIQS_RUNTIME_ERROR << " Overflow of counter ";

  foreach(data.dets[a_level], [this, s](std::pair<time_pt, int> const& x, std::pair<time_pt, int> const& y, double M) {
   // beta-periodicity is implicit in the argument, just fix the sign properly
   this->g_tau[closest_mesh_pt(double(y.first - x.first))](y.second, x.second) += (y.first >= x.first ? real(s) : -real(s)) * M;
  });
 }
 // ---------------------------------------------

 void collect_results(boost::mpi::communicator const& c) {

  boost::mpi::all_reduce(c, z, z, std::c14::plus<>());
  boost::mpi::all_reduce(c, num, num, std::c14::plus<>());
  average_sign = z / num;
  // Need a copy, because all_reduce wants default-constructible types
  auto g_tau_copy = make_clone(g_tau);
  boost::mpi::all_reduce(c, g_tau_copy, g_tau_copy, std::c14::plus<>());
  g_tau = g_tau_copy / (-real(z) * data.config.beta() * g_tau_copy.mesh().delta());
 }
};
}
