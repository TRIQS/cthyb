/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2018, The Simons Foundation & H. U.R. Strand
 * Author: H. U.R. Strand
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

#include "./G2_iw_acc.hpp"

namespace triqs_cthyb {

  namespace G2_iw {

    template <G2_channel Channel>
    measure_G2_iw_base<Channel>::measure_G2_iw_base(std::optional<G2_iw_t> &G2_iw_opt,
                                                       qmc_data const &data,
                                                       G2_measures_t const &G2_measures)
       : data(data), average_sign(0), G2_measures(G2_measures) {

      const double beta = data.config.beta();

      order           = G2_measures.params.measure_G2_block_order;
      int n_bosonic   = G2_measures.params.measure_G2_n_bosonic;
      int n_fermionic = G2_measures.params.measure_G2_n_fermionic;

      // Allocate the two-particle Green's function
      {
        mesh::imfreq mesh_f{beta, Fermion, n_fermionic};
        mesh::imfreq mesh_b{beta, Boson, n_bosonic};

        mesh::prod<imfreq, imfreq, imfreq> mesh_fff{mesh_f, mesh_f, mesh_f};
        mesh::prod<imfreq, imfreq, imfreq> mesh_bff{mesh_b, mesh_f, mesh_f};

        if (Channel == G2_channel::AllFermionic)
          G2_iw_opt = make_block2_gf(mesh_fff, G2_measures.gf_struct, order);
        else
          G2_iw_opt = make_block2_gf(mesh_bff, G2_measures.gf_struct, order);

        G2_iw.rebind(*G2_iw_opt);
        G2_iw() = 0;
      }

      // Allocate temporary two-frequency matrix M
      {
        if (Channel == G2_channel::AllFermionic) { // Smaller mesh possible in AllFermionic
          mesh::imfreq iw_mesh_large{beta, Fermion, 3 * n_fermionic};
          mesh::imfreq iw_mesh_small{beta, Fermion, n_fermionic};
          M_mesh = mesh::prod<imfreq, imfreq>{iw_mesh_large, iw_mesh_small};
        } else {
          int nfreq = n_bosonic + n_fermionic;
          mesh::imfreq iw_mesh{beta, Fermion, nfreq};
          M_mesh = mesh::prod<imfreq, imfreq>{iw_mesh, iw_mesh};
        }

        // Initialize intermediate scattering matrix
        M = block_gf{M_mesh, G2_measures.gf_struct};
      }
    }

    template <G2_channel Channel>
    void accumulate_impl_AABB(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij,
                              M_t const &M_kl);
    template <G2_channel Channel>
    void accumulate_impl_ABBA(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij,
                              M_t const &M_kl);

    template <G2_channel Channel> void measure_G2_iw_base<Channel>::accumulate_G2(mc_weight_t s) {

      s *= data.atomic_reweighting;
      average_sign += s;
      
      timer_G2.start();
      for (auto &m : G2_measures()) {
        auto G2_iw_block = G2_iw(m.b1.index(), m.b2.idx);
        bool diag_block  = (m.b1.index() == m.b2.idx);
        if (order == block_order::AABB || diag_block)
          accumulate_impl_AABB<Channel>(G2_iw_block, s, M(m.b1.index()), M(m.b2.idx));
        if (order == block_order::ABBA || diag_block)
          accumulate_impl_ABBA<Channel>(G2_iw_block, s, M(m.b1.index()), M(m.b2.idx));
      }
      timer_G2.stop();
    }

    template <G2_channel Channel>
    void measure_G2_iw_base<Channel>::collect_results(mpi::communicator const &com) {

      average_sign = mpi::all_reduce(average_sign, com);
      G2_iw        = mpi::all_reduce(G2_iw, com);

      G2_iw = G2_iw / (real(average_sign) * data.config.beta());

      if (com.rank() == 0) {
        std::cout << "measure/G2_iw: timer_M  = " << double(timer_M) << "\n";
        std::cout << "measure/G2_iw: timer_G2 = " << double(timer_G2) << "\n";
      }
    }

    // -- Particle-hole

    template <>
    void accumulate_impl_AABB<G2_channel::PH>(G2_iw_t::g_t::view_type G2, mc_weight_t s,
                                              M_t const &M_ij, M_t const &M_kl) {

      //G2(w, n1, n2)(i, j, k, l) << G2(w, n1, n2)(i, j, k, l) + s * M_ij(n1, n1 + w)(i, j) * M_kl(n2 + w, n2)(k, l);

      for (const auto &[w, n1, n2] : G2.mesh())
        for (const auto [i, j, k, l] : G2.target_indices())
          G2[w, n1, n2](i, j, k, l) += s * M_ij[n1, n1 + w](i, j) * M_kl[n2 + w, n2](k, l);
    }

    template <>
    void accumulate_impl_ABBA<G2_channel::PH>(G2_iw_t::g_t::view_type G2, mc_weight_t s,
                                              M_t const &M_il, M_t const &M_kj) {

      //G2(w, n1, n2)(i, j, k, l) << G2(w, n1, n2)(i, j, k, l) - s * M_il(n1, n2)(i, l) * M_kj(n2 + w, n1 + w)(k, j);

      for (const auto &[w, n1, n2] : G2.mesh())
        for (const auto [i, j, k, l] : G2.target_indices())
          G2[w, n1, n2](i, j, k, l) -= s * M_il[n1, n2](i, l) * M_kj[n2 + w, n1 + w](k, j);
    }

    // -- Particle-particle

    template <>
    void accumulate_impl_AABB<G2_channel::PP>(G2_iw_t::g_t::view_type G2, mc_weight_t s,
                                              M_t const &M_ij, M_t const &M_kl) {

      //G2(w, n1, n2)(i, j, k, l) << G2(w, n1, n2)(i, j, k, l) + s * M_ij(n1, w - n2)(i, j) * M_kl(w - n1, n2)(k, l);

      for (const auto &[w, n1, n2] : G2.mesh())
        for (const auto [i, j, k, l] : G2.target_indices())
          G2[w, n1, n2](i, j, k, l) += s * M_ij[n1, w - n2](i, j) * M_kl[w - n1, n2](k, l);
    }

    template <>
    void accumulate_impl_ABBA<G2_channel::PP>(G2_iw_t::g_t::view_type G2, mc_weight_t s,
                                              M_t const &M_il, M_t const &M_kj) {

      //G2(w, n1, n2)(i, j, k, l) << G2(w, n1, n2)(i, j, k, l) - s * M_il(n1, n2)(i, l) * M_kj(w - n1, w - n2)(k, j);

      for (const auto &[w, n1, n2] : G2.mesh())
        for (const auto [i, j, k, l] : G2.target_indices())
          G2[w, n1, n2](i, j, k, l) -= s * M_il[n1, n2](i, l) * M_kj[w - n1, w - n2](k, j);
    }

    // -- Fermionic

    template <>
    void accumulate_impl_AABB<G2_channel::AllFermionic>(G2_iw_t::g_t::view_type G2, mc_weight_t s,
                                                        M_t const &M_ij, M_t const &M_kl) {

      //G2(n1, n2, n3)(i, j, k, l) << G2(n1, n2, n3)(i, j, k, l) + s * M_ij(n2, n1)(j, i) * M_kl(n1 + n3 - n2, n3)(l, k);

      const auto &iw_mesh = std::get<0>(G2.mesh());
      using mesh_point_t  = typename std::remove_reference<decltype(iw_mesh)>::type::mesh_point_t;

      for (const auto &[n1, n2, n3] : G2.mesh()) {
        mesh_point_t n4{iw_mesh, n1.index() + n3.idx - n2.idx};
        for (const auto [i, j, k, l] : G2.target_indices())
          G2[n1, n2, n3](i, j, k, l) += s * M_ij[n2, n1](j, i) * M_kl[n4, n3](l, k);
      }
    }

    template <>
    void accumulate_impl_ABBA<G2_channel::AllFermionic>(G2_iw_t::g_t::view_type G2, mc_weight_t s,
                                                        M_t const &M_il, M_t const &M_kj) {

      //G2(n1, n2, n3)(i, j, k, l) << G2(n1, n2, n3)(i, j, k, l) - s * M_il(n1 + n3 - n2, n1)(l, i) * M_kj(n2, n3)(j, k);

      const auto &iw_mesh = std::get<0>(G2.mesh());
      using mesh_point_t  = typename std::remove_reference<decltype(iw_mesh)>::type::mesh_point_t;

      for (const auto &[n1, n2, n3] : G2.mesh()) {
        mesh_point_t n4{iw_mesh, n1.index() + n3.idx - n2.idx};
        for (const auto [i, j, k, l] : G2.target_indices())
          G2[n1, n2, n3](i, j, k, l) -= s * M_il[n4, n1](l, i) * M_kj[n2, n3](j, k);
      }
    }

    template class measure_G2_iw_base<G2_channel::AllFermionic>;
    template class measure_G2_iw_base<G2_channel::PP>;
    template class measure_G2_iw_base<G2_channel::PH>;
    
  } // namespace G2_iw
} // namespace triqs_cthyb
