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
#include "nfft_buf.hpp"

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

 auto do_transform = [this](det_type const& det, details::nfft_matrix_t & nfft_matrix) {
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
  nfft_matrix_ab.reset();
  do_transform(det_A, nfft_matrix_ab);
  nfft_matrix_ab.transform();

  nfft_matrix_cd.reset();
  do_transform(det_B, nfft_matrix_cd);
  nfft_matrix_cd.transform();

  if(Channel == PH)
   g2(iw_,inu_,inup_)(a_,b_,c_,d_) << g2(iw_,inu_,inup_)(a_,b_,c_,d_) +
    s * nfft_matrix_ab(a_,b_)(-inu_,inu_+iw_) * nfft_matrix_cd(c_,d_)(-inup_-iw_,inup_);
  else
   g2(iw_,inu_,inup_)(a_,b_,c_,d_) << g2(iw_,inu_,inup_)(a_,b_,c_,d_) +
    s * nfft_matrix_ab(a_,b_)(-inu_,iw_-inup_) * nfft_matrix_cd(c_,d_)(-iw_+inu_,inup_);
 }
 if(Order == ABBA || diag_block) {
  nfft_matrix_ad.reset();
  do_transform(det_A, nfft_matrix_ad);
  nfft_matrix_ad.transform();

  nfft_matrix_cb.reset();
  do_transform(det_B, nfft_matrix_cb);
  nfft_matrix_cb.transform();

  if(Channel == PH)
   g2(iw_,inu_,inup_)(a_,b_,c_,d_) << g2(iw_,inu_,inup_)(a_,b_,c_,d_) -
    s * nfft_matrix_ad(a_,d_)(-inu_,inup_) * nfft_matrix_cb(c_,b_)(-inup_-iw_,inu_+iw_);
  else
   g2(iw_,inu_,inup_)(a_,b_,c_,d_) << g2(iw_,inu_,inup_)(a_,b_,c_,d_) -
    s * nfft_matrix_ad(a_,d_)(-inu_,inup_) * nfft_matrix_cb(c_,b_)(-iw_+inu_,iw_-inup_);
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

nfft_matrix_t::nfft_matrix_t(int size1, int size2, double beta, int n_freq1, int n_freq2) :
 data(size1, size2), beta(beta), max_n_tau(0) {
 data() = data_t{{}, res_gf_t{{{beta, Fermion, n_freq1}, {beta, Fermion, n_freq2}}, {}}};
}

void nfft_matrix_t::reset() {
 foreach(data, [this](int i1, int i2){
  auto & d = data(i1,i2);
  d.input.clear();
 });
 max_n_tau = 0;
}

void nfft_matrix_t::push_back(std::pair<time_pt,int> const& x, std::pair<time_pt,int> const& y, dcomplex fxy) {
 auto & input = data(x.second, y.second).input;
 input.emplace_back(std::array<double,2>{double(x.first), double(y.first)}, fxy);
 max_n_tau = std::max(max_n_tau, input.size());
}

nfft_matrix_t::res_gf_view nfft_matrix_t::operator()(int n1, int n2) const {
 return data(n1,n2).result;
}

void nfft_matrix_t::transform() {
 nfft_buf_t<2> buf(data(0,0).result.data(), max_n_tau, beta, do_checks);

 foreach(data, [this,&buf](int i1, int i2) {
  auto & d = data(i1,i2);
  if(d.input.empty()) {
   d.result() = 0;
   return;
  }
  buf.set_fiw_arr(d.result.data());
  for(auto const& x : d.input) buf.push_back(x.first, x.second);
  if(!buf.is_empty()) buf.flush();
 });
}

}

}
