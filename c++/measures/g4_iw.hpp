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

#include <vector>
#include <triqs/mpi/base.hpp>
#include <triqs/statistics/histograms.hpp>

#include "../nfft_array.hpp"
#include "../qmc_data.hpp"

#include "util.hpp"

namespace cthyb {

  using namespace triqs::arrays;
  
  // Measure the two-particle Green's function in Matsubara frequency
  template <g4_channel Channel> struct measure_g4_iw {

    using M_block_type = block_gf<cartesian_product<imfreq, imfreq>, matrix_valued>;
    using M_type = M_block_type::g_t;
    
    qmc_data const &data;
    g4_iw_t::view_type g4_iw;
    mc_weight_t average_sign;
    block_order order;
    g4_measures_t g4_measures;

    M_block_type M;
    array<nfft_array_t<2, 2>, 1> M_nfft;

    measure_g4_iw(std::optional<g4_iw_t> & g4_iw_opt, qmc_data const &data, g4_measures_t const & g4_measures);
    void accumulate(mc_weight_t s);
    void collect_results(triqs::mpi::communicator const &c);

    inline void accumulate_impl_AABB(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const & M_ab, M_type const & M_cd);
    inline void accumulate_impl_ABBA(g4_iw_t::g_t::view_type g4, mc_weight_t s, M_type const & M_ad, M_type const & M_cb);

  };

}
