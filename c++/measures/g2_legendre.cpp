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
#include "./g2.hpp"

#ifndef NDEBUG
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#endif

namespace cthyb {

  template <g2_channel Channel, block_order Order>
  measure_g2_legendre<Channel, Order>::measure_g2_legendre(int A, int B, g2_view_type g2, qmc_data const &data)
     : A(A), B(B), diag_block(A == B), g2(g2), data(data), size_A(data.delta[A].target_shape()[0]), size_B(data.delta[B].target_shape()[0]) {
    // TODO
  }

  template <g2_channel Channel, block_order Order> void measure_g2_legendre<Channel, Order>::accumulate(mc_weight_t s) {
    // TODO
  }

  template <g2_channel Channel, block_order Order> void measure_g2_legendre<Channel, Order>::collect_results(triqs::mpi::communicator const &c) {
    // TODO
  }

  template class measure_g2_legendre<PP, AABB>;
  template class measure_g2_legendre<PP, ABBA>;
  template class measure_g2_legendre<PH, AABB>;
  template class measure_g2_legendre<PH, ABBA>;
}
