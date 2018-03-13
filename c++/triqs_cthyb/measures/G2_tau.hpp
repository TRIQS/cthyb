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

#include <triqs/gfs.hpp>

#include "../qmc_data.hpp"
#include "util.hpp"

namespace triqs_cthyb {

  using namespace triqs::gfs;

  // Measure imaginary time Green's function (all blocks)
  struct measure_G2_tau {

    public:
    measure_G2_tau(std::optional<G2_tau_t> &G2_tau_opt, qmc_data const &data, G2_measures_t const &G2_measures);

    void accumulate(mc_weight_t sign);
    void collect_results(triqs::mpi::communicator const &comm);

    private:
    qmc_data const &data;
    G2_tau_t::view_type G2_tau;
    mc_weight_t average_sign;
    block_order order;
    G2_measures_t G2_measures;
  };

} // namespace triqs_cthyb
