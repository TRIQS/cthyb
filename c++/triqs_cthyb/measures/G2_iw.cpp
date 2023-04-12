/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2018, The Simons Foundation
 * Author: H. U.R. Strand
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

#include <nda/nda.hpp>
#include <triqs/utility/itertools.hpp>

#include "./G2_iw.hpp"

namespace triqs_cthyb {

  using namespace G2_iw;

  template <G2_channel Channel>
  measure_G2_iw<Channel>::measure_G2_iw(std::optional<G2_iw_t> &G2_iw_opt, qmc_data const &data,
                                        G2_measures_t const &G2_measures)
     : measure_G2_iw_base<Channel>(G2_iw_opt, data, G2_measures) {

    // Accumulation buffer for scattering matrix
    for (auto const &m : M) {
      auto norb1       = static_cast<size_t>(m.target_shape()[0]);
      auto norb2       = static_cast<size_t>(m.target_shape()[1]);
      size_t nfreq_pts = static_cast<size_t>(std::get<0>(M_mesh.components()).size());
      M_block_arr.push_back(M_arr_t(norb1, norb2, nfreq_pts, nfreq_pts));
    }
  }

  template <G2_channel Channel> void measure_G2_iw<Channel>::accumulate(mc_weight_t s) {

    if (true)
      accumulate_M_opt(); // FLOPS Optimized scattering matrix accumulation
    else {

      // ---------------------------------------------------------------
      // Base line scattering matrix computation in two frequencies

      auto M_ww_fill = [this](det_type const &det, M_t &M_ww) {
        const double beta = this->data.config.beta();
        foreach (det, [&M_ww, beta](op_t const &x, op_t const &y, det_scalar_t M_xy) {
          // insert accumulation
          double t1 = double(x.first);
          double t2 = double(y.first);
          for (auto [w1, w2] : M_ww.mesh()) { M_ww[w1, w2](x.second, y.second) += exp((beta - t1) * w1) * M_xy * std::exp(t2 * w2); }
        })
          ;
      };

      timer_M.start();
      // Intermediate M matrices for all blocks
      M() = 0;
      for (auto bidx : range(M.size())) { M_ww_fill(data.dets[bidx], M[bidx]); }
      timer_M.stop();

    } // end accumulate_M

    // ---------------------------------------------------------------
    // Recombine products of scattering matrices
    // to accumulate two particle quantities

    accumulate_G2(s);
  }

  template <G2_channel Channel> void measure_G2_iw<Channel>::accumulate_M_opt() {

    // ---------------------------------------------------------------
    // Scattering matrix accumulation with minimal exponent evaluations
    // using product relations for imaginary time + frequenc exponents

    const double beta    = data.config.beta();
    const double pi_beta = M_PI / beta;

    auto M_arr_fill = [pi_beta, beta](det_type const &det, M_arr_t &M_arr, M_mesh_t const &M_mesh) {
      foreach (det,
               [&M_mesh, &M_arr, pi_beta, beta](op_t const &x, op_t const &y, det_scalar_t M_xy) {
                 double t1 = double(x.first);
                 double t2 = double(y.first);

                 const auto &mesh1 = std::get<0>(M_mesh.components());
                 const auto &mesh2 = std::get<1>(M_mesh.components());

                 std::complex<double> dWt1(0., 2 * pi_beta * (beta - t1));
                 std::complex<double> dWt2(0., 2 * pi_beta * t2);

                 auto dexp1 = std::exp(dWt1);
                 auto dexp2 = std::exp(dWt2);

                 auto exp1 = std::exp(dWt1 * (mesh1.first_index() + 0.5));

                 for (auto const i1 : range(M_arr.shape()[2])) {

                   auto exp2      = std::exp(dWt2 * (mesh2.first_index() + 0.5));
                   auto exp1_M_xy = exp1 * M_xy;

                   for (auto const i2 : range(M_arr.shape()[3])) {

                     M_arr(x.second, y.second, i1, i2) += exp1_M_xy * exp2;
                     exp2 *= dexp2;
                   }
                   exp1 *= dexp1;
                 }
               })
        ;
    };

    timer_M.start();

    // Intermediate M matrices for all blocks
    for (auto bidx : range(M_block_arr.size())) {
      M_block_arr[bidx]() = 0;
      M_arr_fill(data.dets[bidx], M_block_arr[bidx], M_mesh);
    }

    // Reshuffle the accumulated scattering matrix into a Green's function object
    for (auto bidx : range(M_block_arr.size()))
      for (auto [n1, n2, i, j] : product_range(M[bidx].data().shape()))
        M[bidx].data()(n1, n2, i, j) = M_block_arr[bidx](i, j, n1, n2);

    timer_M.stop();
  }

  template class measure_G2_iw<G2_channel::AllFermionic>;
  template class measure_G2_iw<G2_channel::PP>;
  template class measure_G2_iw<G2_channel::PH>;

} // namespace triqs_cthyb
