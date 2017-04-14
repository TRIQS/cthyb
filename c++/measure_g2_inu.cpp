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
#include "measure_g2.hpp"

namespace cthyb {

template<g2_channel Channel, block_order Order>
measure_g2_inu<Channel,Order>::measure_g2_inu(int b1, int b2, g2_iw_inu_inup_block_view g2, qmc_data const& data) :
 b1(b1), b2(b2), g2(g2), data(data) {
// TODO
}

template<g2_channel Channel, block_order Order>
void measure_g2_inu<Channel,Order>::accumulate(mc_weight_t s) {
// TODO
}

template<g2_channel Channel, block_order Order>
void measure_g2_inu<Channel,Order>::collect_results(triqs::mpi::communicator const& c) {
// TODO
}

template class measure_g2_inu<PP,AABB>;
template class measure_g2_inu<PP,ABBA>;
template class measure_g2_inu<PH,AABB>;
template class measure_g2_inu<PH,ABBA>;
}
