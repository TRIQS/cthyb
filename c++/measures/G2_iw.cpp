/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016, P. Seth, I. Krivenko, H. U.R. Strand, M. Ferrero and O. Parcollet
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

#include "./G2_iw.hpp"

namespace cthyb {

  template <G2_channel Channel>
  measure_G2_iw<Channel>::measure_G2_iw(std::optional<G2_iw_t> &G2_iw_opt, qmc_data const &data, G2_measures_t const &G2_measures)
     : data(data), average_sign(0), G2_measures(G2_measures) {

    const double beta = data.config.beta();

    order           = G2_measures.params.measure_G2_block_order;
    int n_bosonic   = G2_measures.params.measure_G2_n_bosonic;
    int n_fermionic = G2_measures.params.measure_G2_n_fermionic;

    // Allocate the two-particle Green's function
    {
      gf_mesh<imfreq> mesh_f{beta, Fermion, n_fermionic};
      gf_mesh<imfreq> mesh_b{beta, Boson, n_bosonic};

      gf_mesh<cartesian_product<imfreq, imfreq, imfreq>> mesh_fff{mesh_f, mesh_f, mesh_f};
      gf_mesh<cartesian_product<imfreq, imfreq, imfreq>> mesh_bff{mesh_b, mesh_f, mesh_f};

      if (Channel == AllFermionic)
        G2_iw_opt = make_block2_gf(mesh_fff, G2_measures.gf_struct, order);
      else
        G2_iw_opt = make_block2_gf(mesh_bff, G2_measures.gf_struct, order);

      G2_iw.rebind(*G2_iw_opt);
      G2_iw() = 0;
    }

    // Allocate temporary NFFT two-frequency matrix M
    {
      int resize_factor = 3; // How much bigger should the large mesh bee???
      int nfreq         = resize_factor * std::max(n_fermionic, n_bosonic - 1 + n_fermionic);
      gf_mesh<imfreq> iw_mesh_large{beta, Fermion, nfreq};
      gf_mesh<cartesian_product<imfreq, imfreq>> M_mesh{iw_mesh_large, iw_mesh_large};

      // Initialize intermediate scattering matrix
      M = make_block_gf(M_mesh, G2_measures.gf_struct);
    }

    // Initialize the nfft_buffers mirroring the matrix M
    {
      M_nfft.resize(M.size());

      for (auto bidx : range(M.size())) {
        std::string bname = M.block_names()[bidx];

        int buf_size = 10; // default
        if (G2_measures.params.nfft_buf_sizes.count(bname)) { buf_size = G2_measures.params.nfft_buf_sizes.at(bname); }
        array<int, 2> buf_sizes{M(bidx).target_shape()};
        buf_sizes() = buf_size;

        M_nfft(bidx) = nfft_array_t<2, 2>(M(bidx).mesh(), M(bidx).data(), buf_sizes);
      }
    }
  }

  template <G2_channel Channel> void measure_G2_iw<Channel>::accumulate(mc_weight_t s) {

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

    for (auto &m : G2_measures()) {
      auto G2_iw_block = G2_iw(m.b1.idx, m.b2.idx);
      bool diag_block  = (m.b1.idx == m.b2.idx);
      if (order == AABB || diag_block) accumulate_impl_AABB(G2_iw_block, s, M(m.b1.idx), M(m.b2.idx));
      if (order == ABBA || diag_block) accumulate_impl_ABBA(G2_iw_block, s, M(m.b1.idx), M(m.b2.idx));
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

  template <> inline void measure_G2_iw<PH>::accumulate_impl_AABB(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_type const &M_ij, M_type const &M_kl) {
    G2(w, n1, n2)(i, j, k, l) << G2(w, n1, n2)(i, j, k, l) //
      + s * M_ij(n1, n1 + w)(i, j) * M_kl(n2 + w, n2)(k, l); // sign in lhs in fft
  }

  template <> inline void measure_G2_iw<PH>::accumulate_impl_ABBA(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_type const &M_il, M_type const &M_kj) {
    G2(w, n1, n2)(i, j, k, l) << G2(w, n1, n2)(i, j, k, l) //
      - s * M_il(n1, n2)(i, l) * M_kj(n2 + w, n1 + w)(k, j); // sign in lhs in fft
  }

  // -- Particle-particle

  template <> inline void measure_G2_iw<PP>::accumulate_impl_AABB(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_type const &M_ij, M_type const &M_kl) {
    G2(w, n1, n2)(i, j, k, l) << G2(w, n1, n2)(i, j, k, l) //
      + s * M_ij(n1, w - n2)(i, j) * M_kl(w - n1, n2)(k, l); // sign in lhs in fft
  }

  template <> inline void measure_G2_iw<PP>::accumulate_impl_ABBA(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_type const &M_il, M_type const &M_kj) {
    G2(w, n1, n2)(i, j, k, l) << G2(w, n1, n2)(i, j, k, l) //
      - s * M_il(n1, n2)(i, l) * M_kj(w - n1, w - n2)(k, j); // sign in lhs in fft
  }

  // -- Fermionic

  template <>
  inline void measure_G2_iw<AllFermionic>::accumulate_impl_AABB(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_type const &M_ij, M_type const &M_kl) {

      int size_ij = M_ij.target_shape()[0];
      int size_kl = M_kl.target_shape()[0];

      for (auto const &n1 : std::get<0>(G2.mesh()))
        for (auto const &n2 : std::get<1>(G2.mesh()))
          for (auto const &n3 : std::get<2>(G2.mesh())) {
            auto mesh = std::get<0>(G2.mesh());
            typename decltype(mesh)::mesh_point_t n4{mesh, n1.index() + n3.index() - n2.index()};
            for (int i : range(size_ij))
              for (int j : range(size_ij))
                for (int k : range(size_kl))
                  for (int l : range(size_kl)) G2[{n1, n2, n3}](i, j, k, l) += s * M_ij[{n2, n1}](j, i) * M_kl[{n4, n3}](l, k);
          }

      //G2(n1, n2, n3)(i, j, k, l) << G2(n1, n2, n3)(i, j, k, l) + s * M_ij(n2, n1)(j, i) * M_kl(n1 + n3 - n2, n3)(l, k);
  }

  template <>
  inline void measure_G2_iw<AllFermionic>::accumulate_impl_ABBA(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_type const &M_il, M_type const &M_kj) {

    int size_il = M_il.target_shape()[0];
    int size_kj = M_kj.target_shape()[0];

    for (auto const &n1 : std::get<0>(G2.mesh()))
      for (auto const &n2 : std::get<1>(G2.mesh()))
        for (auto const &n3 : std::get<2>(G2.mesh())) {
          auto mesh = std::get<0>(G2.mesh());
          typename decltype(mesh)::mesh_point_t n4{mesh, n1.index() + n3.index() - n2.index()};
          for (int i : range(size_il))
            for (int j : range(size_kj))
              for (int k : range(size_kj))
                for (int l : range(size_il)) G2[{n1, n2, n3}](i, j, k, l) += -s * M_il[{n4, n1}](l, i) * M_kj[{n2, n3}](j, k);
        }

    //G2(n1, n2, n3)(i, j, k, l) << G2(n1, n2, n3)(i, j, k, l) - s * M_il(n1 + n3 - n2, n1)(l, i) * M_kj(n2, n3)(j, k);
  }

  // --

  template <G2_channel Channel> void measure_G2_iw<Channel>::collect_results(triqs::mpi::communicator const &com) {
    average_sign = mpi_all_reduce(average_sign, com);
    G2_iw        = mpi_all_reduce(G2_iw, com);
    for (auto const &G2_iw_block : G2_iw) G2_iw_block /= real(average_sign) * data.config.beta();
  }

  template class measure_G2_iw<AllFermionic>;
  template class measure_G2_iw<PP>;
  template class measure_G2_iw<PH>;

} // namespace cthyb
