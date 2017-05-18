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

namespace cthyb {

  using namespace triqs::gfs;

  // Measure imaginary time Green's function (all blocks)
  struct measure_g4_tau {

    qmc_data const &data;
    g4_tau_t::view_type g4_tau;
    mc_weight_t average_sign;
    g4_measures_t g4_measures;

    measure_g4_tau(std::optional<g4_tau_t> & g4_tau_opt, qmc_data const &data, g4_measures_t const & g4_measures);

    void accumulate(mc_weight_t sign);
    void collect_results(triqs::mpi::communicator const &comm);
  };
}
