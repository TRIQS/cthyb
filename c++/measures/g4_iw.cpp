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

#include "./g4_iw.hpp"

namespace cthyb {

  template <g4_channel Channel>
  measure_g4_iw<Channel>::measure_g4_iw(std::optional<g4_iw_t> &g4_iw_opt, qmc_data const &data, g4_measures_t const &g4_measures)
     : data(data), average_sign(0), g4_measures(g4_measures) {

    const double beta = data.config.beta();

    order           = g4_measures.params.measure_g4_block_order;
    int n_bosonic   = g4_measures.params.measure_g4_n_bosonic;
    int n_fermionic = g4_measures.params.measure_g4_n_fermionic;

    // Allocate the two-particle Green's function
    {
      gf_mesh<imfreq> mesh_f{beta, Fermion, n_fermionic};
      gf_mesh<imfreq> mesh_b{beta, Boson, n_bosonic};

      gf_mesh<cartesian_product<imfreq, imfreq, imfreq>> mesh_fff{mesh_f, mesh_f, mesh_f};
      gf_mesh<cartesian_product<imfreq, imfreq, imfreq>> mesh_bff{mesh_b, mesh_f, mesh_f};

      if (Channel == AllFermionic)
        g4_iw_opt = make_block2_gf(mesh_fff, g4_measures.gf_struct);
      else
        g4_iw_opt = make_block2_gf(mesh_bff, g4_measures.gf_struct);

      g4_iw.rebind(*g4_iw_opt);
      g4_iw() = 0;
    }

    // Allocate temporary NFFT two-frequency matrix M
    {
      int resize_factor = 3; // How much bigger should the large mesh bee???
      int nfreq         = resize_factor * std::max(n_fermionic, n_bosonic - 1 + n_fermionic);
      gf_mesh<imfreq> iw_mesh_large{beta, Fermion, nfreq};
      gf_mesh<cartesian_product<imfreq, imfreq>> M_mesh{iw_mesh_large, iw_mesh_large};

      // Initialize intermediate scattering matrix
      M = make_block_gf(M_mesh, g4_measures.gf_struct);
    }

    // Initialize the nfft_buffers mirroring the matrix M
    {
      M_nfft.resize(M.size());

      for (auto bidx : range(M.size())) {
        std::string bname = M.block_names()[bidx];

        int buf_size = 10; // default
        if (g4_measures.params.nfft_buf_sizes.count(bname)) { buf_size = g4_measures.params.nfft_buf_sizes.at(bname); }
        array<int, 2> buf_sizes{M(bidx).target_shape()};
        buf_sizes() = buf_size;

        M_nfft(bidx) = nfft_array_t<2, 2>(M(bidx).mesh(), M(bidx).data(), buf_sizes);
      }
    }
  }

  template <g4_channel Channel> void measure_g4_iw<Channel>::accumulate(mc_weight_t s) {

    s *= data.atomic_reweighting;
    average_sign += s;

    auto nfft_fill = [this](det_type const &det, nfft_array_t<2, 2> &nfft_matrix) {
      foreach (det, [&nfft_matrix](op_t const &x, op_t const &y, det_scalar_t M) {
        nfft_matrix.push_back({double(x.first), double(y.first)}, {x.second, y.second}, M);
      })
        ;
    };

    // Intermediate M matrices for all blocks
    M() = 0;
    for (auto bidx : range(M.size())) {
      nfft_fill(data.dets[bidx], M_nfft(bidx));
      M_nfft(bidx).flush();
    }

    for (auto &m : g4_measures()) {
      auto g4_iw_block = g4_iw(m.b1.idx, m.b2.idx);
      bool diag_block  = (m.b1.idx == m.b2.idx);
      if (order == AABB || diag_block) accumulate_impl_AABB(g4_iw_block, s, M(m.b1.idx), M(m.b2.idx));
      if (order == ABBA || diag_block) accumulate_impl_AABB(g4_iw_block, s, M(m.b1.idx), M(m.b2.idx));
    }
  }

  // Frequency placeholders
  clef::placeholder<0> iw_;
  clef::placeholder<1> inu_;
  clef::placeholder<2> inup_;

  clef::placeholder<0> w1;
  clef::placeholder<1> w2;
  clef::placeholder<2> w3;

  // Index placeholders
  clef::placeholder<3> a_;
  clef::placeholder<4> b_;
  clef::placeholder<5> c_;
  clef::placeholder<6> d_;

  template <> inline void measure_g4_iw<PH>::accumulate_impl_AABB(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const &M_ab, M_type const &M_cd) {
    g4(iw_, inu_, inup_)(a_, b_, c_, d_) << g4(iw_, inu_, inup_)(a_, b_, c_, d_)
          + s * M_ab(-inu_, inu_ + iw_)(a_, b_) * M_cd(-inup_ - iw_, inup_)(c_, d_);
  }

  template <> inline void measure_g4_iw<PP>::accumulate_impl_AABB(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const &M_ab, M_type const &M_cd) {
    g4(iw_, inu_, inup_)(a_, b_, c_, d_) << g4(iw_, inu_, inup_)(a_, b_, c_, d_)
          + s * M_ab(-inu_, iw_ - inup_)(a_, b_) * M_cd(-iw_ + inu_, inup_)(c_, d_);
  }

  template <>
  inline void measure_g4_iw<AllFermionic>::accumulate_impl_AABB(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const &M_ab, M_type const &M_cd) {
    g4(w1, w2, w3)(a_, b_, c_, d_) << g4(w1, w2, w3)(a_, b_, c_, d_) + s * M_ab(w2, w1)(b_, a_) * M_cd(w1 + w3 - w2, w3)(d_, c_);
  }

  // -------------------------

  template <> inline void measure_g4_iw<PH>::accumulate_impl_ABBA(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const &M_ad, M_type const &M_cb) {

    g4(iw_, inu_, inup_)(a_, b_, c_, d_) << g4(iw_, inu_, inup_)(a_, b_, c_, d_)
          - s * M_ad(-inu_, inup_)(a_, d_) * M_cb(-inup_ - iw_, inu_ + iw_)(c_, b_);
  }

  template <> inline void measure_g4_iw<PP>::accumulate_impl_ABBA(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const &M_ad, M_type const &M_cb) {

    g4(iw_, inu_, inup_)(a_, b_, c_, d_) << g4(iw_, inu_, inup_)(a_, b_, c_, d_)
          - s * M_ad(-inu_, inup_)(a_, d_) * M_cb(-iw_ + inu_, iw_ - inup_)(c_, b_);
  }

  template <>
  inline void measure_g4_iw<AllFermionic>::accumulate_impl_ABBA(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const &M_ad, M_type const &M_cb) {

    g4(w1, w2, w3)(a_, b_, c_, d_) << g4(w1, w2, w3)(a_, b_, c_, d_) - s * M_ad(w1 + w3 - w2, w1)(d_, a_) * M_cb(w2, w3)(b_, c_);
  }

  template <g4_channel Channel> void measure_g4_iw<Channel>::collect_results(triqs::mpi::communicator const &c) {
    average_sign = mpi_all_reduce(average_sign, c);

    g4_iw = mpi_all_reduce(g4_iw, c);

    for (auto const &g4_iw_block : g4_iw) g4_iw_block /= real(average_sign) * data.config.beta();
  }

  template class measure_g4_iw<AllFermionic>;
  template class measure_g4_iw<PP>;
  template class measure_g4_iw<PH>;

} // namespace cthyb
