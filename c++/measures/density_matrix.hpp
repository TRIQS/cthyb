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
#include "./qmc_data.hpp"

namespace cthyb {

  struct measure_density_matrix {
    qmc_data const &data;
    std::vector<matrix_t> &block_dm; // density matrix of each block
    mc_weight_t z = 0;

    measure_density_matrix(qmc_data const &data, std::vector<matrix_t> &density_matrix);
    void accumulate(mc_weight_t s);
    void collect_results(triqs::mpi::communicator const &c);
  };
}
