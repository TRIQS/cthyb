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

namespace triqs_cthyb {

  struct measure_average_sign {

    qmc_data const &data;
    mc_weight_t &average_sign;
    mc_weight_t sign, z;

    measure_average_sign(qmc_data const &_data, mc_weight_t &_average_sign) : data(_data), average_sign(_average_sign) {
      average_sign = 1.0;
      z            = 0;
      sign         = 0;
    }
    // --------------------

    void accumulate(mc_weight_t s) {

      sign += s * data.atomic_reweighting;
      z += std::abs(data.atomic_reweighting);
    }
    // ---------------------------------------------

    void collect_results(mpi::communicator const &c) {

      z            = mpi::all_reduce(z, c);
      sign         = mpi::all_reduce(sign, c);
      average_sign = sign / z;
    }
  };
} // namespace triqs_cthyb
