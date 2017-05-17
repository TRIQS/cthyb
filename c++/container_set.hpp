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
    g_tau_t::view_type G_tau() { return *g_tau; }

    // Intermediate object to accumulate g(tau), either real or complex
    std::optional<g_tau_g_target_t> g_tau_accum;
    
    /// Legendre Green's function
    std::optional<g_l_t> g_l;
    g_l_t::view_type G_l() { return *g_l; }

    /// Two-particle Green's function (three fermionic imaginary times)
    std::optional<g4_tau_t> g4_tau;
    g4_tau_t::view_type G2_tau() { return *g4_tau; }

    // Two-particle Green's function (three fermionic matsubaras)
    std::optional<g4_iw_t> g4_iw;
    g4_iw_t::view_type G2_inu() { return *g4_iw;}

    // Two-particle Green's function (three fermionic matsubaras)
    std::optional<g4_iw_t> g4_iw_ph;
    g4_iw_t::view_type G2_iw_inu_inup_ph() { return *g4_iw_ph;}

    // Two-particle Green's function (three fermionic matsubaras)
    std::optional<g4_iw_t> g4_iw_pp;
    g4_iw_t::view_type G2_iw_inu_inup_pp() { return *g4_iw_pp;}
    
    /*
    /// Write containers to hdf5 file
    friend void h5_write(triqs::h5::group h5group, std::string subgroup_name, container_set_t const &c) {

      triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

      h5_write(grp, "G_tau", c.g_tau);
      h5_write(grp, "G_l", c.g_l);
    }

    /// Read containers from hdf5 file
    friend void h5_read(triqs::h5::group h5group, std::string subgroup_name, container_set_t &c) {
      triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

      h5_read(grp, "G_tau", c.g_tau);
      h5_read(grp, "G_l", c.g_l);
    }
    */
    
  }; // struct container_set_t

} // namespace cthyb
