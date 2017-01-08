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
#pragma once

#include <triqs/mpi/base.hpp>
#include "qmc_data.hpp"

namespace cthyb {

enum g2_channel {PP, PH};

// Measure one block of G^2(iomega,inu,inu')
template<g2_channel Channel, block_order Order> struct measure_g2_inu {

 qmc_data const& data;
 g2_iw_inu_inup_block_view g2;
 int b1; // Index of block A within gf_struct
 int b2; // Index of block B within gf_struct

 measure_g2_inu(int b1, int b2, g2_iw_inu_inup_block_view g2, qmc_data const& data);
 void accumulate(mc_weight_t s);
 void collect_results(triqs::mpi::communicator const& c);

};

// Measure one block of G^2(iomega,l,l')
template<g2_channel Channel, block_order Order> struct measure_g2_legendre {

 qmc_data const& data;
 g2_iw_l_lp_block_view g2;
 int b1; // Index of block A within gf_struct
 int b2; // Index of block B within gf_struct

 measure_g2_legendre(int b1, int b2, g2_iw_l_lp_block_view g2, qmc_data const& data);
 void accumulate(mc_weight_t s);
 void collect_results(triqs::mpi::communicator const& c);

};

}
