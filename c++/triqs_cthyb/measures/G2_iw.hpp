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

#include "G2_iw_acc.hpp"

namespace triqs_cthyb {

  // Measure the two-particle Green's function in Matsubara frequency
  template <G2_channel Channel> class measure_G2_iw : public G2_iw::measure_G2_iw_base<Channel> {

    public:
    measure_G2_iw(std::optional<G2_iw_t> &G2_iw_opt, qmc_data const &data,
                  G2_measures_t const &G2_measures);
    void accumulate(mc_weight_t s);
    void accumulate_M_opt(mc_weight_t s);

    using B = G2_iw::measure_G2_iw_base<Channel>;
    using B::collect_results;
    
    private:
    G2_iw::M_block_arr_t M_block_arr;
    using B::M, B::M_mesh, B::G2_measures, B::data, B::timer_M, B::accumulate_G2;
  };

} // namespace triqs_cthyb
