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
        g4_iw_opt = make_block2_gf(mesh_fff, g4_measures.gf_struct, order);
      else
        g4_iw_opt = make_block2_gf(mesh_bff, g4_measures.gf_struct, order);

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
      const double beta = this->data.config.beta();
      foreach (det, [&nfft_matrix, beta](op_t const &x, op_t const &y, det_scalar_t M) {
        nfft_matrix.push_back({beta - double(x.first), double(y.first)}, {x.second, y.second}, M);
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
      if (order == ABBA || diag_block) accumulate_impl_ABBA(g4_iw_block, s, M(m.b1.idx), M(m.b2.idx));
    }
  }

  // Index placeholders
  clef::placeholder<0> i;
  clef::placeholder<1> j;
  clef::placeholder<2> k;
  clef::placeholder<3> l;

  // Frequency placeholders
  clef::placeholder<4> w;
  clef::placeholder<5> n1;
  clef::placeholder<6> n2;
  clef::placeholder<7> n3;

  // -- Particle-hole

  template <> inline void measure_g4_iw<PH>::accumulate_impl_AABB(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const &M_ij, M_type const &M_kl) {
    g4(w, n1, n2)(i, j, k, l) << g4(w, n1, n2)(i, j, k, l) + s * M_ij(-n1, n1 + w)(i, j) * M_kl(-n2 - w, n2)(k, l);
  }

  template <> inline void measure_g4_iw<PH>::accumulate_impl_ABBA(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const &M_il, M_type const &M_kj) {
    g4(w, n1, n2)(i, j, k, l) << g4(w, n1, n2)(i, j, k, l) - s * M_il(-n1, n2)(i, l) * M_kj(-n2 - w, n1 + w)(k, j);
  }

  // -- Particle-particle

  template <> inline void measure_g4_iw<PP>::accumulate_impl_AABB(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const &M_ij, M_type const &M_kl) {
    g4(w, n1, n2)(i, j, k, l) << g4(w, n1, n2)(i, j, k, l) + s * M_ij(-n2, w - n2)(i, j) * M_kl(-w + n1, n2)(k, l);
  }

  template <> inline void measure_g4_iw<PP>::accumulate_impl_ABBA(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const &M_il, M_type const &M_kj) {
    g4(w, n1, n2)(i, j, k, l) << g4(w, n1, n2)(i, j, k, l) - s * M_il(-n1, n2)(i, l) * M_kj(-w + n2, w - n2)(k, j);
  }

  // -- Fermionic

  template <>
  inline void measure_g4_iw<AllFermionic>::accumulate_impl_AABB(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const &M_ij, M_type const &M_kl) {

      int size_ij = M_ij.target_shape()[0];
      int size_kl = M_kl.target_shape()[0];

      for (auto const &n1 : std::get<0>(g4.mesh()))
        for (auto const &n2 : std::get<1>(g4.mesh()))
          for (auto const &n3 : std::get<2>(g4.mesh())) {
            auto mesh = std::get<0>(g4.mesh());
            typename decltype(mesh)::mesh_point_t n4{mesh, n1.index() + n3.index() - n2.index()};
            for (int i : range(size_ij))
              for (int j : range(size_ij))
                for (int k : range(size_kl))
                  for (int l : range(size_kl)) g4[{n1, n2, n3}](i, j, k, l) += s * M_ij[{n2, n1}](j, i) * M_kl[{n4, n3}](l, k);
          }

      //g4(n1, n2, n3)(i, j, k, l) << g4(n1, n2, n3)(i, j, k, l) + s * M_ij(n2, n1)(j, i) * M_kl(n1 + n3 - n2, n3)(l, k);
  }

  template <>
  inline void measure_g4_iw<AllFermionic>::accumulate_impl_ABBA(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const &M_il, M_type const &M_kj) {

    int size_il = M_il.target_shape()[0];
    int size_kj = M_kj.target_shape()[0];

    for (auto const &n1 : std::get<0>(g4.mesh()))
      for (auto const &n2 : std::get<1>(g4.mesh()))
        for (auto const &n3 : std::get<2>(g4.mesh())) {
          auto mesh = std::get<0>(g4.mesh());
          typename decltype(mesh)::mesh_point_t n4{mesh, n1.index() + n3.index() - n2.index()};
          for (int i : range(size_il))
            for (int j : range(size_kj))
              for (int k : range(size_kj))
                for (int l : range(size_il)) g4[{n1, n2, n3}](i, j, k, l) += -s * M_il[{n4, n1}](l, i) * M_kj[{n2, n3}](j, k);
        }

    //g4(n1, n2, n3)(i, j, k, l) << g4(n1, n2, n3)(i, j, k, l) - s * M_il(n1 + n3 - n2, n1)(l, i) * M_kj(n2, n3)(j, k);
  }

  // --

  template <g4_channel Channel> void measure_g4_iw<Channel>::collect_results(triqs::mpi::communicator const &com) {
    average_sign = mpi_all_reduce(average_sign, com);
    g4_iw        = mpi_all_reduce(g4_iw, com);
    for (auto const &g4_iw_block : g4_iw) g4_iw_block /= real(average_sign) * data.config.beta();
  }

  template class measure_g4_iw<AllFermionic>;
  template class measure_g4_iw<PP>;
  template class measure_g4_iw<PH>;

} // namespace cthyb
