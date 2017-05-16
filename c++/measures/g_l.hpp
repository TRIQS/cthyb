/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, H. U.R. Strand, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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

  // Measure Legendre Green's function (all blocks)
  struct measure_g_l {

    qmc_data const &data;
    g_l_view_t g_l;
    mc_weight_t average_sign;

    measure_g_l(g_l_opt_t & g_l_opt, qmc_data const &data, int n_l, gf_struct_t gf_struct) : data(data), average_sign(0) {
      // Allocate storage in g_l_opt
      g_l_opt = make_block_gf(g_l_t::g_t::mesh_t{data.config.beta(), Fermion, static_cast<size_t>(n_l)}, gf_struct);
      g_l.rebind(*g_l_opt);
      g_l() = 0.0;
    }

    // --------------------

    void accumulate(mc_weight_t s) {
      s *= data.atomic_reweighting;
      average_sign += s;

      double beta = data.config.beta();
      auto Tn = triqs::utility::legendre_generator();

      for (auto block_idx : range(g_l.size())) {

        foreach (data.dets[block_idx], [this, s, block_idx, beta, &Tn](op_t const &x, op_t const &y, det_scalar_t M) {

          double poly_arg = 2 * double(y.first - x.first) / beta - 1.0;
          Tn.reset(poly_arg);

          auto val = (y.first >= x.first ? s : -s) * M;

          for (auto l : g_l[block_idx].mesh()) {
            // Evaluate all polynomial orders
            this->g_l[block_idx][l](y.second, x.second) += val * Tn.next();
          }
        })
          ;
      } // for block_idx
    }
    // ---------------------------------------------

    void collect_results(triqs::mpi::communicator const &c) {

      average_sign = mpi_all_reduce(average_sign, c);
      g_l = mpi_all_reduce(g_l, c);

      double beta = data.config.beta();

      for (auto block_idx : range(g_l.size())) {
        for (auto l : g_l[block_idx].mesh()) {
          /// Normalize polynomial coefficients with basis overlap
          g_l[block_idx][l] *= -(sqrt(2.0 * l + 1.0) / (real(average_sign) * beta));
        }
	matrix<double> id(g_l[block_idx].target_shape());
	id() = 1.0; // this creates an unit matrix
        enforce_discontinuity(g_l[block_idx], id);
      }
    }
    // ---------------------------------------------
  };
}
