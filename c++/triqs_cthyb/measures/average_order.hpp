/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2021, Simons Foundation
 *    author: N. Wentzell
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

namespace triqs_cthyb {

  /// Measure of the average perturbation order
  struct measure_average_order {

    measure_average_order(qmc_data const &_data, double &_average_order) : data(_data), average_order(_average_order) { average_order = 0.0; }

    void accumulate(mc_weight_t) {
      average_order += data.config.size() / 2;
      ++N;
    }

    void collect_results(mpi::communicator const &comm) {
      N = mpi::all_reduce(N, comm);

      // Reduce and normalize
      average_order = mpi::all_reduce(average_order, comm);
      average_order = average_order / N;
    }

    private:
    // The Monte-Carlo configuration
    qmc_data const &data;

    // Reference to double for accumulation
    double &average_order;

    // Accumulation counter
    long long N = 0;
  };

} // namespace triqs_cthyb
