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

#include "types.hpp"

namespace triqs_cthyb {

  // Containers for measurements
  struct container_set_t {

    // -- Single particle Green's functions

    /// Single-particle Green's function :math:`G(\tau)` in imaginary time.
    std::optional<G_tau_t> G_tau;

    /// Intermediate Green's function to accumulate g(tau), either real or complex
    std::optional<G_tau_G_target_t> G_tau_accum;

    /// Violation of the fundamental Green function property G(tau)[i,j] = G(tau)*[j,i] after the measurement
    std::optional<G_tau_G_target_t> asymmetry_G_tau;

    /// Single-particle Green's function :math:`G(\tau)` in imaginary time.
    std::optional<G_iw_t> G_iw_direct;

    /// Single-particle Green's function :math:`G_l` in Legendre polynomial representation.
    std::optional<G_l_t> G_l;

    /// General operator Green's function :math:`O(\tau)` in imaginary time.
    std::optional<gf<imtime, scalar_valued>> O_tau;

    // -- Two-particle Green's functions

    /// Two-particle Green's function :math:`G^{(2)}(\tau_1,\tau_2,\tau_3)` (three Fermionic imaginary times)
    std::optional<G2_tau_t> G2_tau;

    /// Two-particle Green's function :math:`G^{(2)}(i\nu,i\nu',i\nu'')` (three Fermionic frequencies)
    std::optional<G2_iw_t> G2_iw;

    /// Two-particle Green's function :math:`G^{(2)}(i\nu,i\nu',i\nu'')` (three Fermionic frequencies)
    std::optional<G2_iw_t> G2_iw_nfft;

    /// Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the pp-channel (one bosonic matsubara and two fermionic)
    std::optional<G2_iw_t> G2_iw_pp;

    /// Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the pp-channel (one bosonic matsubara and two fermionic)
    std::optional<G2_iw_t> G2_iw_pp_nfft;

    /// Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the ph-channel (one bosonic matsubara and two fermionic)
    std::optional<G2_iw_t> G2_iw_ph;

    /// Two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the ph-channel (one bosonic matsubara and two fermionic)
    std::optional<G2_iw_t> G2_iw_ph_nfft;

    /// Two-particle Green's function :math:`G^{(2)}(i\omega,l,l')` in the pp-channel (one bosonic matsubara and two legendre)
    std::optional<G2_iwll_t> G2_iwll_pp;

    /// Two-particle Green's function :math:`G^{(2)}(i\omega,l,l')` in the ph-channel (one bosonic matsubara and two legendre)
    std::optional<G2_iwll_t> G2_iwll_ph;

    /// Histogram of the total perturbation order
    std::optional<histogram> perturbation_order_total;

    /// Histograms of the perturbation order for each block
    std::optional<histo_map_t> perturbation_order;

    /// Function that writes all containers to hdf5 file
    friend void h5_write(h5::group h5group, std::string subgroup_name, container_set_t const &c);

    /// Function that reads all containers to hdf5 file
    friend void h5_read(h5::group h5group, std::string subgroup_name, container_set_t &c);

  }; // struct container_set_t

} // namespace triqs_cthyb
