/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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
  struct measure_g {

    qmc_data const &data;
    gf_view<imtime, g_target_t> g_tau;
    int a_level;
    mc_weight_t z;
    int64_t num;
    mc_weight_t average_sign;

    measure_g(int a_level, gf_view<imtime, g_target_t> g_tau, qmc_data const &data);
    void accumulate(mc_weight_t s);
    void collect_results(triqs::mpi::communicator const &c);
  };
}
