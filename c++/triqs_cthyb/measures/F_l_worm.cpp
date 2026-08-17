/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2026
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

#include "./F_l_worm.hpp"

namespace triqs_cthyb {

  using namespace triqs::gfs;
  using namespace triqs::mesh;

  namespace {

    int worm_time_order_sign(qmc_data const &data) {
      auto const &worm = data.worm;
      int exponent     = 0;

      for (auto const &[tau, op] : data.config) {
        auto in_interval = (worm.tau_Q >= worm.tau_cdag) ? (tau > worm.tau_cdag && tau <= worm.tau_Q)
                                                         : (tau > worm.tau_cdag || tau <= worm.tau_Q);
        if (in_interval) ++exponent;
      }

      return (exponent % 2 == 0 ? 1 : -1);
    }

  } // namespace

  measure_F_l_worm::measure_F_l_worm(std::optional<G_l_t> &F_l_worm_opt, qmc_data const &data, int n_l,
                                     gf_struct_t const &gf_struct, double worm_eta)
     : data(data), worm_eta(worm_eta), sign_Z(0) {
    F_l_worm_opt = block_gf<legendre>{{data.config.beta(), Fermion, n_l}, gf_struct};
    F_l_worm.rebind(*F_l_worm_opt);
    F_l_worm() = 0.0;
  }

  void measure_F_l_worm::accumulate(mc_weight_t s) {
    s *= data.atomic_reweighting;

    if (data.worm.in_Z()) {
      sign_Z += s;
      return;
    }

    auto const &worm = data.worm;
    double beta      = data.config.beta();
    double poly_arg  = 2 * double(worm.tau_Q - worm.tau_cdag) / beta - 1.0;
    auto Tn          = triqs::utility::legendre_generator();
    Tn.reset(poly_arg);

    auto val = worm_time_order_sign(data) * s;
    for (auto l : F_l_worm[worm.block].mesh()) { F_l_worm[worm.block][l](worm.inner_Q, worm.inner_cdag) += val * Tn.next(); }
  }

  void measure_F_l_worm::collect_results(mpi::communicator const &c) {
    F_l_worm = mpi::all_reduce(F_l_worm, c);
    sign_Z   = mpi::all_reduce(sign_Z, c);

    double beta = data.config.beta();
    for (auto &F_l_block : F_l_worm)
      for (auto l : F_l_block.mesh()) F_l_block[l] *= sqrt(2.0 * l.index() + 1.0) / (real(sign_Z) * beta * worm_eta);
  }

} // namespace triqs_cthyb
