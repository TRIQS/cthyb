/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, H. U.R. Strand
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

#include <optional>

#include <types.hpp>

namespace cthyb {

  /// Containers for measurements
  struct container_set_t {

    /// Imaginary-time Green's function
    std::optional<g_tau_t> g_tau;

    /// Accumulated :math:`G(\tau)` in imaginary time.
    g_tau_t::view_type G_tau() { return *g_tau; }

    // Intermediate object to accumulate g(tau), either real or complex
    std::optional<g_tau_g_target_t> g_tau_accum;
    
    /// Legendre Green's function
    std::optional<g_l_t> g_l;

    /// Accumulated :math:`G_l` in Legendre polynomials representation.
    g_l_t::view_type G_l() { return *g_l; }

    /// Two-particle Green's function (three fermionic imaginary times)
    std::optional<g4_tau_t> g4_tau;

    /// Accumulated two-particle Green's function :math:`G^{(2)}(\tau_1,\tau_2,\tau_3)`
    g4_tau_t::view_type G2_tau() { return *g4_tau; }

    // Two-particle Green's function (three fermionic matsubaras)
    std::optional<g4_iw_t> g4_iw;

    /// Accumulated two-particle Green's function :math:`G^{(2)}(i\nu,i\nu',i\nu'')`
    g4_iw_t::view_type G2_inu() { return *g4_iw;}

    // Two-particle Green's function (three fermionic matsubaras)
    std::optional<g4_iw_t> g4_iw_pp;

    /// Accumulated two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the pp-channel.
    g4_iw_t::view_type G2_iw_inu_inup_pp() { return *g4_iw_pp;}

    // Two-particle Green's function (three fermionic matsubaras)
    std::optional<g4_iw_t> g4_iw_ph;

    /// Accumulated two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the ph-channel.
    g4_iw_t::view_type G2_iw_inu_inup_ph() { return *g4_iw_ph;}

    // Two-particle Green's function (one bosonic matsubara and two legendre)
    std::optional<g4_wll_t> g4_wll_pp;

    /// Accumulated two-particle Green's function :math:`G^{(2)}(i\omega,l,l')` in the pp-channel.
    g4_wll_t::view_type G2_iw_l_lp_pp() { return *g4_wll_pp;}

    // Two-particle Green's function (one bosonic matsubara and two legendre)
    std::optional<g4_wll_t> g4_wll_ph;

    /// Accumulated two-particle Green's function :math:`G^{(2)}(i\omega,l,l')` in the ph-channel.
    g4_wll_t::view_type G2_iw_l_lp_ph() { return *g4_wll_ph;}
    
  }; // struct container_set_t

} // namespace cthyb
