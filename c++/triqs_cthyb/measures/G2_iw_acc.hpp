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
#pragma once

#include "./G2_iw.hpp"

namespace triqs_cthyb {
  namespace G2_iw {

    template<G2_channel Channel> void accumulate_impl_AABB(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);
    template<G2_channel Channel> void accumulate_impl_ABBA(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);

    /*
    template<G2_channel Channel> void accumulate_impl_AABB_opt(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);
    template<G2_channel Channel> void accumulate_impl_ABBA_opt(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);

    template<G2_channel Channel> void accumulate_impl_AABB_w0(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);
    template<G2_channel Channel> void accumulate_impl_ABBA_w0(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);

    template<G2_channel Channel> void accumulate_impl_AABB_opt_w0(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);
    template<G2_channel Channel> void accumulate_impl_ABBA_opt_w0(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);
    */
    
  } // namespace G2_iw
} // namespace triqs_cthyb
