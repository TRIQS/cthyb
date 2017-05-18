/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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

#include "./g2.hpp"

#ifndef NDEBUG
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#endif

namespace cthyb {

  template <g2_channel Channel>
  measure_g2_legendre<Channel>::measure_g2_legendre(int A, int B, g2_view_type g2, qmc_data const &data, int buf_size_A, int buf_size_B,
                                                    block_order &order)
     : A(A), B(B), g2(g2), data(data), n_l(std::get<1>(g2.mesh()).size()), order(order), average_sign(0) {

    g2() = 0;

    int size_A = data.delta[A].target_shape()[0];
    int size_B = data.delta[B].target_shape()[0];

    int n_iw = (std::get<0>(g2.mesh()).size() + 1) / 2;

    const double beta = data.config.beta();

    int s1 = size_A;
    int s2 = order == AABB ? size_A : size_B;
    int s3 = order == AABB ? size_B : size_A;
    int s4 = size_B;

    array<int, 6> buf_sizes(n_l, n_l, s1, s2, s3, s4);
    //buf_sizes() = size_A * size_B;
    buf_sizes() = 10;

    std::cout << "n_l, n_iw = " << n_l << ", " << n_iw << "\n";
    std::cout << "buf_target_shape = " << buf_sizes.shape() << "\n";
    std::cout << "g2.data().shape() = " << g2.data().shape() << "\n";

    gf_mesh<imfreq> mesh{beta, Boson, n_iw};
    nfft_buf = nfft_array_t<1, 6>{mesh, g2.data(), buf_sizes};
  }

  template <g2_channel Channel> void measure_g2_legendre<Channel>::accumulate(mc_weight_t s) {

    s *= data.atomic_reweighting;
    average_sign += s;

    if (data.dets[A].size() == 0 || data.dets[B].size() == 0) return;
    
    double beta = data.config.beta();
    
    auto accumulate_impl = [&](op_t const &i, op_t const &j, op_t const &k, op_t const &l, det_scalar_t val) {
      tilde_p_gen p_l1_gen(beta), p_l2_gen(beta);
      double dtau = setup_times(p_l1_gen, p_l2_gen, i, j, k, l);
      for (int l1 : range(n_l) ) {
        double p_l1 = p_l1_gen.next();
        for (int l2 : range(n_l) ) {
          double p_l2 = p_l2_gen.next();
          mini_vector<int, 6> vec{l1, l2, i.second, j.second, k.second, l.second};
          nfft_buf.push_back({dtau}, vec, val * p_l1 * p_l2);
        }
      }
    };

    bool diag_block = (A == B);
    if (order == AABB || diag_block) {
      foreach (data.dets[A], [&](op_t const &i, op_t const &j, det_scalar_t M_ij) {
        foreach (data.dets[B], [&](op_t const &k, op_t const &l, det_scalar_t M_kl) {
          accumulate_impl(i, j, k, l, s * M_ij * M_kl); // Accumulate in legendre-nfft buffer
        })
          ;
      })
        ;
    }
    if (order == ABBA || diag_block) {
      foreach (data.dets[A], [&](op_t const &i, op_t const &l, det_scalar_t M_il) {
        foreach (data.dets[B], [&](op_t const &k, op_t const &j, det_scalar_t M_kj) {
          accumulate_impl(i, j, k, l, -s * M_il * M_kj); // Accumulate in legendre-nfft buffer
        })
          ;
      })
        ;
    }
  }

  template <>
  double measure_g2_legendre<PH>::setup_times(tilde_p_gen &p_l1_gen, tilde_p_gen &p_l2_gen, op_t const &i, op_t const &j, op_t const &k,
                                              op_t const &l) {
    p_l1_gen.reset(i.first, j.first);
    p_l2_gen.reset(k.first, l.first);
    double dtau = 0.5 * double(i.first + j.first - k.first - l.first);
    return dtau;
  }

  template <>
  double measure_g2_legendre<PP>::setup_times(tilde_p_gen &p_l1_gen, tilde_p_gen &p_l2_gen, op_t const &i, op_t const &j, op_t const &k,
                                              op_t const &l) {
    p_l1_gen.reset(i.first, k.first);
    p_l2_gen.reset(j.first, l.first);
    double dtau = 0.5 * double(i.first + mult_by_int(j.first, 3) - mult_by_int(k.first, 3) - l.first);
    return dtau;
  }

  template <g2_channel Channel> void measure_g2_legendre<Channel>::collect_results(triqs::mpi::communicator const &c) {

    nfft_buf.flush();
    g2 = mpi_all_reduce(g2, c);

    average_sign = mpi_all_reduce(average_sign, c);

    for (auto l : std::get<1>(g2.mesh().components())) {
      auto _ = var_t{};
      double s = std::sqrt(2 * l + 1);
      g2[_][l][_] *= s;
      g2[_][_][l] *= s * (l % 2 ? 1 : -1);
    }

    g2 = g2 / (real(average_sign) * data.config.beta());
  }

  template class measure_g2_legendre<PP>;
  template class measure_g2_legendre<PH>;
}
