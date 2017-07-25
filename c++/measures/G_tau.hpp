/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, H. U.R. Strand, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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

#include "./qmc_data.hpp"

namespace cthyb {

  using namespace triqs::gfs;

  // Measure imaginary time Green's function (all blocks)
  class measure_G_tau {

    public:
    measure_G_tau(std::optional<G_tau_G_target_t> &G_tau_opt, qmc_data const &data, int n_tau, gf_struct_t gf_struct);
    void accumulate(mc_weight_t s);
    void collect_results(triqs::mpi::communicator const &c);

    private:
    G_tau_G_target_t::view_type G_tau;
    qmc_data const &data;
    mc_weight_t average_sign;
  };

} // namespace cthyb
