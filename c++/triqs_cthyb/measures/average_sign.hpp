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

#include <limits>

namespace triqs_cthyb {

  struct measure_average_sign {

    qmc_data const &data;
    mc_weight_t &average_sign_partition, &average_sign_worm;
    mc_weight_t sign_partition = 0, sign_worm = 0;
    double norm_partition = 0, norm_worm = 0;

    measure_average_sign(qmc_data const &_data, mc_weight_t &_average_sign_partition, mc_weight_t &_average_sign_worm)
       : data(_data), average_sign_partition(_average_sign_partition), average_sign_worm(_average_sign_worm) {
      average_sign_partition = 1.0;
      average_sign_worm      = std::numeric_limits<double>::quiet_NaN();
    }
    // --------------------

    void accumulate(mc_weight_t s) {

      auto weighted_sign = s * data.atomic_reweighting;
      auto norm          = std::abs(data.atomic_reweighting);

      if (data.worm.in_Z()) {
        sign_partition += weighted_sign;
        norm_partition += norm;
      } else {
        sign_worm += weighted_sign;
        norm_worm += norm;
      }
    }
    // ---------------------------------------------

    void collect_results(mpi::communicator const &c) {

      sign_partition = mpi::all_reduce(sign_partition, c);
      sign_worm      = mpi::all_reduce(sign_worm, c);
      norm_partition = mpi::all_reduce(norm_partition, c);
      norm_worm      = mpi::all_reduce(norm_worm, c);

      auto nan               = std::numeric_limits<double>::quiet_NaN();
      average_sign_partition = norm_partition > 0 ? sign_partition / norm_partition : mc_weight_t{nan};
      average_sign_worm      = norm_worm > 0 ? sign_worm / norm_worm : mc_weight_t{nan};
    }
  };
} // namespace triqs_cthyb
