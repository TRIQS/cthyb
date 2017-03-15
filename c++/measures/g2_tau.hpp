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
#pragma once
#include "../qmc_data.hpp"
#include <triqs/gfs.hpp>

namespace cthyb {

  using namespace triqs::gfs;

  // Measure imaginary time Green's function (one block)
  struct measure_g2_tau {

    using g2_view_type = gf_view<cartesian_product<imtime, imtime, imtime>, tensor_valued<4>>;

    qmc_data const &data;
    g2_view_type g2;
    const int A;           // Index of block A within gf_struct
    const int B;           // Index of block B within gf_struct
    const bool diag_block; // A == B
    mc_weight_t z;
    int64_t num;
    mc_weight_t average_sign;

    measure_g2_tau(int A, int B, g2_view_type g2, qmc_data const &data);

    void accumulate(mc_weight_t s);
    void collect_results(triqs::mpi::communicator const &c);
  };
}
