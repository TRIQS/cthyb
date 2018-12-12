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

    namespace {

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

    } // namespace

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
      
      /*
      triqs::utility::timer timer_clef;
      triqs::utility::timer timer_1loop;
      triqs::utility::timer timer_2loop;
      triqs::utility::timer timer_3loop;
      triqs::utility::timer timer_5loop;
      triqs::utility::timer timer_7loop;

      // --

      timer_clef.start();

      G2(w, n1, n2)
         (i, j, k, l) << G2(w, n1, n2)(i, j, k, l)
            - s * M_il(n1, n2)(i, l) * M_kj(n2 + w, n1 + w)(k, j);

      timer_clef.stop();

      // --

      timer_7loop.start();

      const auto &b_mesh = std::get<0>(G2.mesh());
      const auto &f_mesh = std::get<1>(G2.mesh());

      for (const auto &w : b_mesh)
        for (const auto &n1 : f_mesh)
          for (const auto &n2 : f_mesh)
            for (const auto i : range(G2.target_shape()[0]))
              for (const auto j : range(G2.target_shape()[1]))
                for (const auto k : range(G2.target_shape()[2]))
                  for (const auto l : range(G2.target_shape()[3]))
                    G2[w, n1, n2](i, j, k, l) -=
                       s * M_il[n1, n2](i, l) * M_kj[n2 + w, n1 + w](k, j);

      timer_7loop.stop();

      // --

      timer_5loop.start();

      for (const auto &[w, n1, n2] : G2.mesh())
        for (const auto i : range(G2.target_shape()[0]))
          for (const auto j : range(G2.target_shape()[1]))
            for (const auto k : range(G2.target_shape()[2]))
              for (const auto l : range(G2.target_shape()[3]))
                G2[w, n1, n2](i, j, k, l) -= s * M_il[n1, n2](i, l) * M_kj[n2 + w, n1 + w](k, j);

      timer_5loop.stop();

      // --

      timer_3loop.start();

      for (const auto &w : b_mesh)
        for (const auto &n1 : f_mesh)
          for (const auto &n2 : f_mesh)
            G2[w, n1, n2](0, 0, 0, 0) -= s * M_il[n1, n2](0, 0) * M_kj[n2 + w, n1 + w](0, 0);

      timer_3loop.stop();

      // --

      timer_2loop.start();

      */
      
      for (const auto &[w, n1, n2] : G2.mesh())
        for (const auto [i, j, k, l] : G2.target_indices())
          G2[w, n1, n2](i, j, k, l) -= s * M_il[n1, n2](i, l) * M_kj[n2 + w, n1 + w](k, j);

      /*
      timer_2loop.stop();

      // --

      timer_1loop.start();

      for (const auto &[w, n1, n2] : G2.mesh()) // 2 times slower than 3lvl nested loop!?!?!
        G2[w, n1, n2](0, 0, 0, 0) -= s * M_il[n1, n2](0, 0) * M_kj[n2 + w, n1 + w](0, 0);

      timer_1loop.stop();

      // --

      std::cout << "clef, 7, 5, 3, 2, 1loop, 2loop/7loop, 1loop/3loop, "
                   "2loop/5loop = "
                << double(timer_clef) << ", " << double(timer_7loop) << ", " << double(timer_5loop)
                << ", " << double(timer_3loop) << ", " << double(timer_2loop) << ", "
                << double(timer_1loop) << ", " << double(timer_2loop) / double(timer_7loop) << ", "
                << double(timer_1loop) / double(timer_3loop) << ", "
                << double(timer_2loop) / double(timer_5loop) << ", " << std::endl;
      */
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
        mesh_point_t n4{iw_mesh, n1.index() + n3.index() - n2.index()};
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
        mesh_point_t n4{iw_mesh, n1.index() + n3.index() - n2.index()};
        for (const auto [i, j, k, l] : G2.target_indices())
          G2[n1, n2, n3](i, j, k, l) -= s * M_il[n4, n1](l, i) * M_kj[n2, n3](j, k);
      }
    }

  } // namespace G2_iw
} // namespace triqs_cthyb
