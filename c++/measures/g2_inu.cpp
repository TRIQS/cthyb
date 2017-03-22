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
  measure_g2_inu<Channel, Order>::measure_g2_inu(int A, int B, g2_view_type g2, qmc_data const &data, int buf_size_A, int buf_size_B)
     : A(A), B(B), diag_block(A == B), g2(g2), data(data) {

    int size_A = data.delta[A].target_shape()[0];
    int size_B = data.delta[B].target_shape()[0];

    int n_iw  = (std::get<0>(g2.mesh()).size() + 1) / 2;
    int n_inu = std::get<1>(g2.mesh()).size() / 2;

    g2() = 0;
    z    = 0;
    num  = 0;

    const double beta = data.config.beta();

    array<int, 2> buf_sizes_A(size_A, size_A), buf_sizes_B(size_B, size_B);
    buf_sizes_A() = buf_size_A;
    buf_sizes_B() = buf_size_B;

    if (Order == AABB || diag_block) {
      nfft_matrix_ab = nfft_array_t<2, 2>({{beta, Fermion, n_inu}, {beta, Fermion, n_iw - 1 + n_inu}}, {size_A, size_A}, buf_sizes_A);
      nfft_matrix_cd = nfft_array_t<2, 2>({{beta, Fermion, n_iw - 1 + n_inu}, {beta, Fermion, n_inu}}, {size_B, size_B}, buf_sizes_B);
    }
    if (Order == ABBA || diag_block) {
      nfft_matrix_ad = nfft_array_t<2, 2>({{beta, Fermion, n_inu}, {beta, Fermion, n_inu}}, {size_A, size_A}, buf_sizes_A);
      nfft_matrix_cb = nfft_array_t<2, 2>({{beta, Fermion, n_iw - 1 + n_inu}, {beta, Fermion, n_iw - 1 + n_inu}}, {size_B, size_B}, buf_sizes_B);
    }
  }

  template <g2_channel Channel, block_order Order> void measure_g2_inu<Channel, Order>::accumulate(mc_weight_t s) {
    num += 1;
    if (num < 0) TRIQS_RUNTIME_ERROR << "Overflow of counter";

    s *= data.atomic_reweighting;
    z += s;

    auto const &det_A = data.dets[A];
    auto const &det_B = data.dets[B];
    if (det_A.size() == 0 || det_B.size() == 0) return;

    auto nfft_fill = [this](det_type const &det, nfft_array_t<2, 2> &nfft_matrix) {
      foreach (det,
               [&nfft_matrix](std::pair<time_pt, int> const &x, std::pair<time_pt, int> const &y, det_scalar_t M) {
               nfft_matrix.push_back({double(x.first), double(y.first)}, {x.second, y.second}, M); })
        ;
    };

    // Frequency placeholders
    clef::placeholder<0> iw_;
    clef::placeholder<1> inu_;
    clef::placeholder<2> inup_;

    // Index placeholders
    clef::placeholder<3> a_;
    clef::placeholder<4> b_;
    clef::placeholder<5> c_;
    clef::placeholder<6> d_;

    if (Order == AABB || diag_block) {
      nfft_matrix_ab().data()() = .0;
      nfft_fill(det_A, nfft_matrix_ab);
      nfft_matrix_ab.flush();

      nfft_matrix_cd().data()() = .0;
      nfft_fill(det_B, nfft_matrix_cd);
      nfft_matrix_cd.flush();

      if (Channel == PH)
        g2(iw_, inu_, inup_)(a_, b_, c_, d_) << g2(iw_, inu_, inup_)(a_, b_, c_, d_)
              + s * nfft_matrix_ab()(-inu_, inu_ + iw_)(a_, b_) * nfft_matrix_cd()(-inup_ - iw_, inup_)(c_, d_);
      else
        g2(iw_, inu_, inup_)(a_, b_, c_, d_) << g2(iw_, inu_, inup_)(a_, b_, c_, d_)
              + s * nfft_matrix_ab()(-inu_, iw_ - inup_)(a_, b_) * nfft_matrix_cd()(-iw_ + inu_, inup_)(c_, d_);
    }
    if (Order == ABBA || diag_block) {
      nfft_matrix_ad().data()() = .0;
      nfft_fill(det_A, nfft_matrix_ad);
      nfft_matrix_ad.flush();

      nfft_matrix_cb().data()() = .0;
      nfft_fill(det_B, nfft_matrix_cb);
      nfft_matrix_cb.flush();

      if (Channel == PH)
        g2(iw_, inu_, inup_)(a_, b_, c_, d_) << g2(iw_, inu_, inup_)(a_, b_, c_, d_)
              - s * nfft_matrix_ad()(-inu_, inup_)(a_, d_) * nfft_matrix_cb()(-inup_ - iw_, inu_ + iw_)(c_, b_);
      else
        g2(iw_, inu_, inup_)(a_, b_, c_, d_) << g2(iw_, inu_, inup_)(a_, b_, c_, d_)
              - s * nfft_matrix_ad()(-inu_, inup_)(a_, d_) * nfft_matrix_cb()(-iw_ + inu_, iw_ - inup_)(c_, b_);
    }
  }

  template <g2_channel Channel, block_order Order> void measure_g2_inu<Channel, Order>::collect_results(triqs::mpi::communicator const &c) {
    z  = mpi_all_reduce(z, c);
    g2 = mpi_all_reduce(g2, c);
    g2 = g2 / (real(z) * data.config.beta());
  }

  template class measure_g2_inu<PP, AABB>;
  template class measure_g2_inu<PP, ABBA>;
  template class measure_g2_inu<PH, AABB>;
  template class measure_g2_inu<PH, ABBA>;
}
