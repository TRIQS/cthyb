/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016, P. Seth, I. Krivenko, H. U.R. Strand, M. Ferrero and O. Parcollet
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
#include <triqs/experimental/nfft_array.hpp>

#include "../qmc_data.hpp"

#include "util.hpp"

namespace cthyb {

  using namespace triqs::arrays;
  using namespace triqs::experimental;
  
  // Measure the two-particle Green's function in Matsubara frequency
  template <G2_channel Channel> struct measure_G2_iw {

    using M_block_type = block_gf<cartesian_product<imfreq, imfreq>, matrix_valued>;
    using M_type = M_block_type::g_t;
    
    qmc_data const &data;
    G2_iw_t::view_type G2_iw;
    mc_weight_t average_sign;
    block_order order;
    G2_measures_t G2_measures;

    M_block_type M;
    array<nfft_array_t<2, 2>, 1> M_nfft;

    measure_G2_iw(std::optional<G2_iw_t> & G2_iw_opt, qmc_data const &data, G2_measures_t const & G2_measures);
    void accumulate(mc_weight_t s);
    void collect_results(triqs::mpi::communicator const &c);

    inline void accumulate_impl_AABB(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_type const & M_ab, M_type const & M_cd);
    inline void accumulate_impl_ABBA(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_type const & M_ad, M_type const & M_cb);

  };

}
