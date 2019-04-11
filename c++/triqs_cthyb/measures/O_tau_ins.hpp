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

#include <triqs/gfs.hpp>

#include "../qmc_data.hpp"

namespace triqs_cthyb {

  using namespace triqs::gfs;

  // Measure imaginary time Green's function (all blocks)
  class measure_O_tau_ins {

    public:
    measure_O_tau_ins(std::optional<gf<imtime, scalar_valued>> &O_tau_opt, qmc_data const &data, int n_tau, many_body_op_t const &op1, many_body_op_t const &op2, mc_tools::random_generator &rng);
    void accumulate(mc_weight_t s);
    void collect_results(mpi::communicator const &c);

    private:
    qmc_data const &data;
    mc_weight_t average_sign;
    gf<imtime, scalar_valued>::view_type O_tau;
    many_body_op_t op1, op2;
    op_desc op1_d, op2_d;
    mc_tools::random_generator &rng;
    
  };

} // namespace triqs_cthyb
