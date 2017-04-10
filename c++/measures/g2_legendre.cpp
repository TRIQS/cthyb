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
#include <triqs/utility/legendre.hpp>

#ifndef NDEBUG
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#endif

namespace cthyb {

  // Generates values of \tilde P_l(x(\tau_1-\tau_2))
  struct tilde_p_gen {
    triqs::utility::legendre_generator l_gen;
    double beta;
    double f;
    tilde_p_gen(double beta) : beta(beta) {}
    void reset(time_pt const& tau1, time_pt const& tau2) {
      l_gen.reset(2 * double(tau1 - tau2) / beta - 1);
      f = tau1 > tau2 ? 1 : -1;
    }
    double next() { return f * l_gen.next(); }
  };

  template <g2_channel Channel, block_order Order>
  measure_g2_legendre<Channel, Order>::measure_g2_legendre(int A, int B, g2_view_type g2, qmc_data const &data, int buf_size_A, int buf_size_B)
     : A(A), B(B), diag_block(A == B), g2(g2), data(data), n_l(std::get<1>(g2.mesh()).size()) {

    int size_A = data.delta[A].target_shape()[0];
    int size_B = data.delta[B].target_shape()[0];

    long n_iw  = (std::get<0>(g2.mesh()).size() + 1) / 2;

    g2() = 0;
    z    = 0;
    num  = 0;

    const double beta = data.config.beta();

    int s1 = size_A;
    int s2 = Order == AABB ? size_A : size_B;
    int s3 = Order == AABB ? size_B : size_A;
    int s4 = size_B;

    array<int, 6> buf_sizes(n_l, n_l, s1, s2, s3, s4);
    buf_sizes() = size_A * size_B;

    nfft_abcd = nfft_array_t<1, 6>({{beta, Fermion, n_iw}}, g2.data(), buf_sizes);
  }

  template <g2_channel Channel, block_order Order> void measure_g2_legendre<Channel, Order>::accumulate(mc_weight_t s) {
    num += 1;
    if (num < 0) TRIQS_RUNTIME_ERROR << "Overflow of counter";

    s *= data.atomic_reweighting;
    z += s;

    auto const &det_A = data.dets[A];
    auto const &det_B = data.dets[B];
    if (det_A.size() == 0 || det_B.size() == 0) return;

    using det_arg_t = std::pair<time_pt, int>;

    auto l_range = range(n_l);
    tilde_p_gen p_l1_gen(data.config.beta()), p_l2_gen(data.config.beta());

    if (Order == AABB || diag_block) {
      foreach (data.dets[A], [&](det_arg_t const &i, det_arg_t const &j, det_scalar_t M_ij) {
        foreach (data.dets[B], [&](det_arg_t const &k, det_arg_t const &l, det_scalar_t M_kl) {
          double dtau;
          if(Channel == PH) {
            p_l1_gen.reset(i.first, j.first);
            p_l2_gen.reset(k.first, l.first);
            dtau = 0.5 * double(i.first + j.first - k.first - l.first);
          } else {
            p_l1_gen.reset(i.first, k.first);
            p_l2_gen.reset(j.first, l.first);
            dtau = 0.5 * double(i.first + mult_by_int(j.first, 3) - mult_by_int(k.first, 3) - l.first);
          }

          for(int l1 : l_range) {
            double p_l1 = p_l1_gen.next();
            for(int l2 : l_range) {
              double p_l2 = p_l2_gen.next();
              nfft_abcd.push_back({dtau},
                                  {l1, l2, i.second, j.second, k.second, l.second},
                                  s * p_l1 * p_l2 * M_ij * M_kl);
            }
          }
        });
      });
    }
    if (Order == ABBA || diag_block) {
      foreach (data.dets[A], [&](det_arg_t const &i, det_arg_t const &l, det_scalar_t M_il) {
        foreach (data.dets[B], [&](det_arg_t const &k, det_arg_t const &j, det_scalar_t M_kj) {
          double dtau;
          if(Channel == PH) {
            p_l1_gen.reset(i.first, j.first);
            p_l2_gen.reset(k.first, l.first);
            dtau = 0.5 * double(i.first + j.first - k.first - l.first);
          } else {
            p_l1_gen.reset(i.first, k.first);
            p_l2_gen.reset(j.first, l.first);
            dtau = 0.5 * double(i.first + mult_by_int(j.first, 3) - mult_by_int(k.first, 3) - l.first);
          }

          for(int l1 : l_range) {
            double p_l1 = p_l1_gen.next();
            for(int l2 : l_range) {
              double p_l2 = p_l2_gen.next();
              nfft_abcd.push_back({dtau},
                                  {l1, l2, i.second, j.second, k.second, l.second},
                                  -s * p_l1 * p_l2 * M_il * M_kj);
            }
          }
        });
      });
    }
  }

  template <g2_channel Channel, block_order Order> void measure_g2_legendre<Channel, Order>::collect_results(triqs::mpi::communicator const &c) {
    nfft_abcd.flush();

    z  = mpi_all_reduce(z, c);
    g2 = mpi_all_reduce(g2, c);

    for(auto l : std::get<1>(g2.mesh().components())) {
      double s = std::sqrt(2 * l + 1);
      g2[var_t()][l][var_t()] *= s;
      g2[var_t()][var_t()][l] *= s * (l % 2 ? 1 : -1);
    }

    g2 = g2 / (real(z) * data.config.beta());
  }

  template class measure_g2_legendre<PP, AABB>;
  template class measure_g2_legendre<PP, ABBA>;
  template class measure_g2_legendre<PH, AABB>;
  template class measure_g2_legendre<PH, ABBA>;
}
