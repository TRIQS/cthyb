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

#include <triqs/gfs.hpp>
#include <triqs/gfs/types.hpp>

#include "./g2.hpp"

#ifndef NDEBUG
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#endif

namespace cthyb {

  template <g2_channel Channel>
  measure_g2_inu<Channel>::measure_g2_inu(int A, int B, g2_view_type g2, qmc_data const &data, int buf_size_A, int buf_size_B, block_order order, gf_struct_t const & gf_struct)
     : A(A), B(B), diag_block(A == B), g2(g2), data(data), average_sign(0), order(order) {

    g2() = 0;

    int n_iw  = (std::get<0>(g2.mesh()).size() + 1) / 2;
    int n_inu = std::get<1>(g2.mesh()).size() / 2;

    int size_A = data.delta[A].target_shape()[0];
    int size_B = data.delta[B].target_shape()[0];

    array<int, 2> buf_sizes_A(size_A, size_A), buf_sizes_B(size_B, size_B);
    buf_sizes_A() = buf_size_A;
    buf_sizes_B() = buf_size_B;

    const double beta = data.config.beta();
    gf_mesh<imfreq> mesh1{beta, Fermion, n_inu};
    gf_mesh<imfreq> mesh2{beta, Fermion, n_iw - 1 + n_inu};
    
    if (order == AABB || diag_block) {

      M_ab = M_type{{mesh1, mesh2}, {size_A, size_A}};
      M_cd = M_type{{mesh2, mesh1}, {size_B, size_B}};

      nfft_M_ab = nfft_array_t<2, 2>(M_ab.mesh(), M_ab.data(), buf_sizes_A);
      nfft_M_cd = nfft_array_t<2, 2>(M_cd.mesh(), M_cd.data(), buf_sizes_B);
    }
    if (order == ABBA || diag_block) {

      M_ad = M_type{{mesh1, mesh1}, {size_A, size_A}};
      M_cb = M_type{{mesh2, mesh2}, {size_B, size_B}};

      nfft_M_ad = nfft_array_t<2, 2>(M_ad.mesh(), M_ad.data(), buf_sizes_A);
      nfft_M_cb = nfft_array_t<2, 2>(M_cb.mesh(), M_cb.data(), buf_sizes_B);
    }

    // -----------------------------------------------

    // Construct Matsubara mesh for temporary Matrix
    int resize_factor = 3; // How much bigger should the large mesh bee???
    int nfreq = resize_factor * std::max(n_iw, n_iw - 1 + n_inu);
    gf_mesh<imfreq> iw_mesh_large{beta, Fermion, nfreq};
    gf_mesh<cartesian_product<imfreq, imfreq>> M_mesh{iw_mesh_large, iw_mesh_large};

    // Initialize intermediate scattering matrix
    M = make_block_gf(M_mesh, gf_struct);
    M_nfft.resize(M.size());
    
    for( auto bidx : range(M.size()) ) {
      // the setting up of buf_sizes is a bit ugly... FIXME
      // use M(bidx).target_shape() directly (yielding a mini_vector<unsigned long, 2>)
      auto target_shape = M(bidx).target_shape();
      array<int, 2> buf_sizes{target_shape[0], target_shape[1]};

      buf_sizes() = 100; // make this an external parameter!
      //buf_sizes() = buf_size;
      
      M_nfft(bidx) = nfft_array_t<2, 2>(M(bidx).mesh(), M(bidx).data(), buf_sizes);
    }
    
  }

  template <g2_channel Channel> void measure_g2_inu<Channel>::accumulate(mc_weight_t s) {

    s *= data.atomic_reweighting;
    average_sign += s;

    auto const &det_A = data.dets[A];
    auto const &det_B = data.dets[B];

    if (det_A.size() == 0 || det_B.size() == 0) return;

    auto nfft_fill = [this](det_type const &det, nfft_array_t<2, 2> &nfft_matrix) {
      foreach (det, [&nfft_matrix](op_t const &x, op_t const &y, det_scalar_t M) {
        nfft_matrix.push_back({double(x.first), double(y.first)}, {x.second, y.second}, M);
      })
        ;
    };

    // Intermediate M matrices for all blocks
    M() = 0;
    for( auto bidx : range(M.size()) ) {
      nfft_fill(data.dets[bidx], M_nfft(bidx));
      M_nfft(bidx).flush();
    }
    
    if (order == AABB || diag_block) {

      /*
      M_ab() = 0;
      nfft_fill(det_A, nfft_M_ab);
      nfft_M_ab.flush();
      
      M_cd() = 0;
      nfft_fill(det_B, nfft_M_cd);
      nfft_M_cd.flush();

      accumulate_impl_AABB(s, M_ab, M_cd);
      */
      accumulate_impl_AABB(s, M(A), M(B));
    }
    if (order == ABBA || diag_block) {
      /*
      M_ad() = 0;
      nfft_fill(det_A, nfft_M_ad);
      nfft_M_ad.flush();

      M_cb() = 0;
      nfft_fill(det_B, nfft_M_cb);
      nfft_M_cb.flush();

      accumulate_impl_ABBA(s, M_ad, M_cb);
      */
      accumulate_impl_AABB(s, M(A), M(B));
    }
  }

  template <> void measure_g2_inu<PH>::accumulate_impl_AABB(mc_weight_t s, M_type const &M_ab, M_type const &M_cd) {
    g2(iw_, inu_, inup_)(a_, b_, c_, d_) << g2(iw_, inu_, inup_)(a_, b_, c_, d_)
          + s * M_ab(-inu_, inu_ + iw_)(a_, b_) * M_cd(-inup_ - iw_, inup_)(c_, d_);
  }

  template <> void measure_g2_inu<PP>::accumulate_impl_AABB(mc_weight_t s, M_type const &M_ab, M_type const &M_cd) {
    g2(iw_, inu_, inup_)(a_, b_, c_, d_) << g2(iw_, inu_, inup_)(a_, b_, c_, d_)
          + s * M_ab(-inu_, iw_ - inup_)(a_, b_) * M_cd(-iw_ + inu_, inup_)(c_, d_);
  }

  template <> void measure_g2_inu<AllFermionic>::accumulate_impl_AABB(mc_weight_t s, M_type const &M_ab, M_type const &M_cd) {
    g2(w1, w2, w3)(a_, b_, c_, d_) << g2(w1, w2, w3)(a_, b_, c_, d_) + s * M_ab(w2, w1)(b_, a_) * M_cd(w1 + w3 - w2, w3)(d_, c_);
  }

  // -------------------------

  template <> inline void measure_g2_inu<PH>::accumulate_impl_ABBA(mc_weight_t s, M_type const &M_ad, M_type const &M_cb) {

    g2(iw_, inu_, inup_)(a_, b_, c_, d_) << g2(iw_, inu_, inup_)(a_, b_, c_, d_)
          - s * M_ad(-inu_, inup_)(a_, d_) * M_cb(-inup_ - iw_, inu_ + iw_)(c_, b_);
  }

  template <> inline void measure_g2_inu<PP>::accumulate_impl_ABBA(mc_weight_t s, M_type const &M_ad, M_type const &M_cb) {

    g2(iw_, inu_, inup_)(a_, b_, c_, d_) << g2(iw_, inu_, inup_)(a_, b_, c_, d_)
          - s * M_ad(-inu_, inup_)(a_, d_) * M_cb(-iw_ + inu_, iw_ - inup_)(c_, b_);
  }

  template <> inline void measure_g2_inu<AllFermionic>::accumulate_impl_ABBA(mc_weight_t s, M_type const &M_ad, M_type const &M_cb) {

    g2(w1, w2, w3)(a_, b_, c_, d_) << g2(w1, w2, w3)(a_, b_, c_, d_) - s * M_ab(w1 + w3 - w2, w1)(d_, a_) * M_cd(w2, w3)(b_, c_);
  }

  template <g2_channel Channel> void measure_g2_inu<Channel>::collect_results(triqs::mpi::communicator const &c) {
    average_sign = mpi_all_reduce(average_sign, c);

    g2 = mpi_all_reduce(g2, c);

    g2 = g2 / (real(average_sign) * data.config.beta());
  }

  template class measure_g2_inu<AllFermionic>;
  template class measure_g2_inu<PP>;
  template class measure_g2_inu<PH>;
}
