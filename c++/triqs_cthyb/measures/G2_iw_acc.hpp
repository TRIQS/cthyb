#pragma once

#include "./G2_iw.hpp"

namespace triqs_cthyb {
  namespace G2_iw {

    template<G2_channel Channel> void accumulate_impl_AABB(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);
    template<G2_channel Channel> void accumulate_impl_ABBA(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);

    template<G2_channel Channel> void accumulate_impl_AABB_opt(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);
    template<G2_channel Channel> void accumulate_impl_ABBA_opt(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);

    template<G2_channel Channel> void accumulate_impl_AABB_w0(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);
    template<G2_channel Channel> void accumulate_impl_ABBA_w0(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);

    template<G2_channel Channel> void accumulate_impl_AABB_opt_w0(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);
    template<G2_channel Channel> void accumulate_impl_ABBA_opt_w0(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl);
    
  } // namespace G2_iw
} // namespace triqs_cthyb
