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

#include "G_l.hpp"

namespace triqs_cthyb {

  using namespace triqs::gfs;
  using namespace triqs::mesh;

  measure_G_l::measure_G_l(std::optional<G_l_t> &G_l_opt, qmc_data const &data, int n_l, gf_struct_t const &gf_struct) : data(data), average_sign(0) {
    G_l_opt = block_gf<legendre>{{data.config.beta(), Fermion, static_cast<size_t>(n_l)}, gf_struct};
    G_l.rebind(*G_l_opt);
    G_l() = 0.0;
  }

  void measure_G_l::accumulate(mc_weight_t s) {
    s *= data.atomic_reweighting;
    average_sign += s;

    double beta = data.config.beta();
    auto Tn     = triqs::utility::legendre_generator();

    for (auto block_idx : range(G_l.size())) {

      foreach (data.dets[block_idx], [this, s, block_idx, beta, &Tn](op_t const &x, op_t const &y, det_scalar_t M) {

        double poly_arg = 2 * double(y.first - x.first) / beta - 1.0;
        Tn.reset(poly_arg);

        auto val = (y.first >= x.first ? s : -s) * M;

        for (auto l : G_l[block_idx].mesh()) {
          // Evaluate all polynomial orders
          this->G_l[block_idx][l](y.second, x.second) += val * Tn.next();
        }
      })
        ;
    } // for block_idx
  }

  void measure_G_l::collect_results(mpi::communicator const &c) {

    average_sign = mpi::all_reduce(average_sign, c);
    G_l          = mpi::all_reduce(G_l, c);

    double beta = data.config.beta();

    for (auto &G_l_block : G_l) {
      for (auto l : G_l_block.mesh()) {
        /// Normalize polynomial coefficients with basis overlap
        G_l_block[l] *= -(sqrt(2.0 * l + 1.0) / (real(average_sign) * beta));
      }
      matrix<double> id(G_l_block.target_shape());
      id() = 1.0; // this creates an unit matrix
      enforce_discontinuity(G_l_block, id);
    }
  }

} // namespace triqs_cthyb
