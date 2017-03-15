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

namespace cthyb {

  using namespace triqs::gfs;

  // Measure Legendre Green's function (one block)
  struct measure_g_legendre {

    qmc_data const &data;
    gf_view<legendre> g_l;
    int a_level;
    double beta;
    mc_weight_t z;
    int64_t num;

    measure_g_legendre(int a_level, gf_view<legendre> g_l, qmc_data const &data) : data(data), g_l(g_l), a_level(a_level), beta(data.config.beta()) {
      g_l() = 0.0;
      z     = 0;
      num   = 0;
    }
    // --------------------

    void accumulate(mc_weight_t s) {
      num += 1;
      if (num < 0) TRIQS_RUNTIME_ERROR << " Overflow of counter ";

      s *= data.atomic_reweighting;
      z += s;

      auto Tn = triqs::utility::legendre_generator();

      foreach (data.dets[a_level], [this, s, &Tn](std::pair<time_pt, int> const &x, std::pair<time_pt, int> const &y, det_scalar_t M) {
        double poly_arg = 2 * double(y.first - x.first) / beta - 1.0;
        Tn.reset(poly_arg);

        auto val = (y.first >= x.first ? s : -s) * M;
        for (auto l : g_l.mesh()) this->g_l[l](y.second, x.second) += val * Tn.next();
      })
        ;
    }
    // ---------------------------------------------

    void collect_results(triqs::mpi::communicator const &c) {

      z                                = mpi_all_reduce(z, c);
      g_l                              = mpi_all_reduce(g_l, c);
      for (auto l : g_l.mesh()) g_l[l] = -(sqrt(2.0 * l + 1.0) / (real(z) * beta)) * g_l[l];

      matrix<double> id(g_l.target_shape());
      id() = 1.0;
      enforce_discontinuity(g_l, id);
    }
  };
}
