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
#include <array>
#include "measure_g2.hpp"

#ifndef NDEBUG
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#endif

namespace cthyb {

template<g2_channel Channel, block_order Order>
measure_g2_inu<Channel,Order>::measure_g2_inu(int A, int B, g2_iw_inu_inup_block_view g2, qmc_data const& data) :
 A(A), B(B), diag_block(A == B), g2(g2), data(data) {

 int size_A = get_target_shape(data.delta[A])[0];
 int size_B = get_target_shape(data.delta[B])[0];

 int n_iw = (std::get<0>(g2.mesh()).size()+1)/2;
 int n_inu = std::get<1>(g2.mesh()).size()/2;

 g2() = 0;
 z = 0;
 num = 0;

 using details::nfft_matrix_t;
 if(Order == AABB || diag_block) {
  nfft_matrix_ab = nfft_matrix_t(size_A, size_A, data.config.beta(), n_inu, n_iw-1+n_inu);
  nfft_matrix_cd = nfft_matrix_t(size_B, size_B, data.config.beta(), n_iw-1+n_inu, n_inu);
 }
 if(Order == ABBA || diag_block) {
  nfft_matrix_ad = nfft_matrix_t(size_A, size_A, data.config.beta(), n_inu, n_inu);
  nfft_matrix_cb = nfft_matrix_t(size_B, size_B, data.config.beta(), n_iw-1+n_inu, n_iw-1+n_inu);
 }
}

template<g2_channel Channel, block_order Order>
void measure_g2_inu<Channel,Order>::accumulate(mc_weight_t s) {
 num += 1;
 if (num < 0) TRIQS_RUNTIME_ERROR << "Overflow of counter";

 s *= data.atomic_reweighting;
 z += s;

 auto const& det_A = data.dets[A];
 auto const& det_B = data.dets[B];
 if(det_A.size() == 0 || det_B.size() == 0) return;

 auto nfft_fill = [this](det_type const& det, details::nfft_matrix_t & nfft_matrix) {
  foreach(det, [&nfft_matrix](std::pair<time_pt,int> const& x, std::pair<time_pt,int> const& y,
                              det_scalar_t M) {
   nfft_matrix.push_back(x,y,M);
  });
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

 if(Order == AABB || diag_block) {
  nfft_matrix_ab.resize_bufs(det_A.size()*det_A.size());
  nfft_fill(det_A, nfft_matrix_ab);
  nfft_matrix_ab.transform();

  nfft_matrix_cd.resize_bufs(det_B.size()*det_B.size());
  nfft_fill(det_B, nfft_matrix_cd);
  nfft_matrix_cd.transform();

  if(Channel == PH)
   g2(iw_,inu_,inup_)(a_,b_,c_,d_) << g2(iw_,inu_,inup_)(a_,b_,c_,d_) +
    s * nfft_matrix_ab()(-inu_,inu_+iw_)(a_,b_) * nfft_matrix_cd()(-inup_-iw_,inup_)(c_,d_);
  else
   g2(iw_,inu_,inup_)(a_,b_,c_,d_) << g2(iw_,inu_,inup_)(a_,b_,c_,d_) +
    s * nfft_matrix_ab()(-inu_,iw_-inup_)(a_,b_) * nfft_matrix_cd()(-iw_+inu_,inup_)(c_,d_);
 }
 if(Order == ABBA || diag_block) {
  nfft_matrix_ad.resize_bufs(det_A.size()*det_A.size());
  nfft_fill(det_A, nfft_matrix_ad);
  nfft_matrix_ad.transform();

  nfft_matrix_cb.resize_bufs(det_B.size()*det_B.size());
  nfft_fill(det_B, nfft_matrix_cb);
  nfft_matrix_cb.transform();

  if(Channel == PH)
   g2(iw_,inu_,inup_)(a_,b_,c_,d_) << g2(iw_,inu_,inup_)(a_,b_,c_,d_) -
    s * nfft_matrix_ad()(-inu_,inup_)(a_,d_) * nfft_matrix_cb()(-inup_-iw_,inu_+iw_)(c_,b_);
  else
   g2(iw_,inu_,inup_)(a_,b_,c_,d_) << g2(iw_,inu_,inup_)(a_,b_,c_,d_) -
    s * nfft_matrix_ad()(-inu_,inup_)(a_,d_) * nfft_matrix_cb()(-iw_+inu_,iw_-inup_)(c_,b_);
 }
}

template<g2_channel Channel, block_order Order>
void measure_g2_inu<Channel,Order>::collect_results(triqs::mpi::communicator const& c) {
 z = mpi_all_reduce(z,c);
 g2 = mpi_all_reduce(g2, c);
 g2 = g2 / (real(z) * data.config.beta());
}

template class measure_g2_inu<PP,AABB>;
template class measure_g2_inu<PP,ABBA>;
template class measure_g2_inu<PH,AABB>;
template class measure_g2_inu<PH,ABBA>;

namespace details {

#ifdef NDEBUG
 bool nfft_matrix_t::do_checks = false;
#else
 bool nfft_matrix_t::do_checks = true;
#endif

nfft_matrix_t::nfft_matrix_t(int size1, int size2, double beta, int n_freq1, int n_freq2) :
 size1(size1), size2(size2), max_n_tau(1),
 result({{{beta, Fermion, n_freq1}, {beta, Fermion, n_freq2}}, {size1, size2}}) {
 buffers.reserve(size1*size2);
 for(int a : range(size1))
 for(int b : range(size2))
  buffers.emplace_back(result.mesh(), max_n_tau, do_checks);
}

void nfft_matrix_t::resize_bufs(int n_tau) {
 if(n_tau > max_n_tau) {
  max_n_tau = n_tau;
  for(auto & buf : buffers) buf = {result.mesh(), max_n_tau, do_checks};
 }
}

void nfft_matrix_t::push_back(std::pair<time_pt,int> const& x, std::pair<time_pt,int> const& y, dcomplex fxy) {
 buffers[x.second*size2 + y.second].push_back({double(x.first), double(y.first)}, fxy);
}

auto nfft_matrix_t::operator()() const -> res_gf_t const& { return result; }

void nfft_matrix_t::transform() {
 for(int a : range(size1))
 for(int b : range(size2)) {
  auto & buf = buffers[a*size2 + b];
  if(buf.is_empty()) {
   slice_target_to_scalar(result, a, b) = 0;
  } else {
   buf.flush();
   buf.fill_gf(slice_target_to_scalar(result, a, b));
  }
 }
}

}

}
