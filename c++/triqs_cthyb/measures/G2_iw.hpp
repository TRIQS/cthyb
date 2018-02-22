/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2018, The Simons Foundation
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

#include <vector>
#include <triqs/mpi/base.hpp>
#include <triqs/utility/timer.hpp> // DEBUG

#include "../qmc_data.hpp"
#include "util.hpp"

namespace triqs_cthyb {

  namespace G2_iw {
    using M_block_t = block_gf<cartesian_product<imfreq, imfreq>, matrix_valued>;
    using M_t       = M_block_t::g_t;
    using M_mesh_t = M_block_t::g_t::mesh_t;

    using M_arr_t = array<std::complex<double>, 4>;
    using M_block_arr_t = std::vector<M_arr_t>;
  }

  using namespace triqs::arrays;
  using namespace G2_iw;

  // Measure the two-particle Green's function in Matsubara frequency
  template <G2_channel Channel> struct measure_G2_iw {

    public:
    measure_G2_iw(std::optional<G2_iw_t> &G2_iw_opt, qmc_data const &data, G2_measures_t const &G2_measures);
    void accumulate(mc_weight_t s);
    void collect_results(triqs::mpi::communicator const &c);

    private:

    qmc_data const &data;
    G2_iw_t::view_type G2_iw;
    mc_weight_t average_sign;
    block_order order;
    G2_measures_t G2_measures;

    M_block_t M;
    M_mesh_t M_mesh;
    M_block_arr_t M_block_arr;

    triqs::utility::timer timer_M;
    triqs::utility::timer timer_G2;    
    
  };

} // namespace triqs_cthyb
