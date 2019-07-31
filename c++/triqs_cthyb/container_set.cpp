/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, H. U.R. Strand, M. Ferrero and O. Parcollet
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

#include "./container_set.hpp"

namespace triqs_cthyb {

  /// Function that writes all containers to hdf5 file
  void h5_write(triqs::h5::group h5group, std::string subgroup_name, container_set_t const &c) {

    triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

    h5_write(grp, "G_tau", c.G_tau);
    h5_write(grp, "G_tau_accum", c.G_tau_accum);
    h5_write(grp, "G_l", c.G_l);
    h5_write(grp, "O_tau", c.O_tau);

    h5_write(grp, "G2_tau", c.G2_tau);
    h5_write(grp, "G2_iw", c.G2_iw);
    h5_write(grp, "G2_iw_nfft", c.G2_iw_nfft);
    h5_write(grp, "G2_iw_pp", c.G2_iw_pp);
    h5_write(grp, "G2_iw_pp_nfft", c.G2_iw_pp_nfft);
    h5_write(grp, "G2_iw_ph", c.G2_iw_ph);
    h5_write(grp, "G2_iw_ph_nfft", c.G2_iw_ph_nfft);
    h5_write(grp, "G2_iwll_pp", c.G2_iwll_pp);
    h5_write(grp, "G2_iwll_ph", c.G2_iwll_ph);
  }

  /// Function that reads all containers to hdf5 file
  void h5_read(triqs::h5::group h5group, std::string subgroup_name, container_set_t &c) {

    triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

    h5_read(grp, "G_tau", c.G_tau);
    h5_read(grp, "G_tau_accum", c.G_tau_accum);
    h5_read(grp, "G_l", c.G_l);
    h5_try_read(grp, "O_tau", c.O_tau);

    h5_read(grp, "G2_tau", c.G2_tau);
    h5_read(grp, "G2_iw", c.G2_iw);
    h5_read(grp, "G2_iw_nfft", c.G2_iw_nfft);
    h5_read(grp, "G2_iw_pp", c.G2_iw_pp);
    h5_read(grp, "G2_iw_pp_nfft", c.G2_iw_pp_nfft);
    h5_read(grp, "G2_iw_ph", c.G2_iw_ph);
    h5_read(grp, "G2_iw_ph_nfft", c.G2_iw_ph_nfft);
    h5_read(grp, "G2_iwll_pp", c.G2_iwll_pp);
    h5_read(grp, "G2_iwll_ph", c.G2_iwll_ph);
  }

} // namespace triqs_cthyb
