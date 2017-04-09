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

  template <g2_channel Channel, block_order Order>
  measure_g2_legendre<Channel, Order>::measure_g2_legendre(int A, int B, g2_view_type g2, qmc_data const &data, int buf_size_A, int buf_size_B)
     : A(A), B(B), diag_block(A == B), g2(g2), data(data) {

    int size_A = data.delta[A].target_shape()[0];
    int size_B = data.delta[B].target_shape()[0];

    int n_iw  = (std::get<0>(g2.mesh()).size() + 1) / 2;
    int n_l   = std::get<1>(g2.mesh()).size();

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

    nfft_tensor_abcd = nfft_array_t<1, 6>({{beta, Fermion, n_iw}}, {n_l, n_l, s1, s2, s3, s4}, buf_sizes);
  }

  template <g2_channel Channel, block_order Order> void measure_g2_legendre<Channel, Order>::accumulate(mc_weight_t s) {
    num += 1;
    if (num < 0) TRIQS_RUNTIME_ERROR << "Overflow of counter";

    s *= data.atomic_reweighting;
    z += s;

    // TODO

    // Frequency/Legendre placeholders
    clef::placeholder<0> iw_;
    clef::placeholder<1> l_;
    clef::placeholder<2> lp_;

    // Index placeholders
    clef::placeholder<3> a_;
    clef::placeholder<4> b_;
    clef::placeholder<5> c_;
    clef::placeholder<6> d_;

    // TODO
    //g2(iw_, l_, lp_)(a_, b_, c_, d_) << g2(iw_, l_, lp_)(a_, b_, c_, d_)
    //                                  + s * nfft_tensor_abcd()(iw_)(l_, lp_, a_, b_, c_, d_);
  }

  template <g2_channel Channel, block_order Order> void measure_g2_legendre<Channel, Order>::collect_results(triqs::mpi::communicator const &c) {
    z  = mpi_all_reduce(z, c);
    g2 = mpi_all_reduce(g2, c);

    for(auto l : std::get<1>(g2.mesh().components())) {
      double s = std::sqrt(2 * l + 1);
      g2[var_t()][l][var_t()] *= s;
      g2[var_t()][var_t()][l] *= s * (l % 2 ? -1 : 1);
    }

    g2 = g2 / (real(z) * data.config.beta());
  }

  template class measure_g2_legendre<PP, AABB>;
  template class measure_g2_legendre<PP, ABBA>;
  template class measure_g2_legendre<PH, AABB>;
  template class measure_g2_legendre<PH, ABBA>;
}
