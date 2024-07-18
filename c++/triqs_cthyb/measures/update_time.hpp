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

  /// Measure of the average update time
  struct measure_update_time {

    measure_update_time(qmc_data &_data, double &_update_time) : data(_data), update_time(_update_time) { update_time = 0.0; }

    void accumulate(mc_weight_t) {
      if (data.updated && !flag) {
        update_time += data.n_acc;
        data.n_acc = 1.;
	++N;
	data.updated = false;
      }
      else
        data.n_acc += 1.;
      if (flag) data.updated = false; // need to set it back to false at the first measurement since it has been set to true during warmup
      flag = false;
    }

    void collect_results(mpi::communicator const &comm) {
      if (N == 0) {
        N = 1;
	update_time = data.n_acc;
      }
	    
      N = mpi::all_reduce(N, comm);

      // Reduce and normalize
      update_time = mpi::all_reduce(update_time, comm);
      update_time = update_time / N;
    }

    private:
    // The Monte-Carlo configuration
    qmc_data &data;

    // Flag to indicate whether or not this is the first measure
    bool flag = true;

    // Reference to double for accumulation
    double &update_time;

    // Accumulation counter
    long long N = 0;
  };

} // namespace triqs_cthyb
