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
#include <triqs/utility/legendre.hpp>
#include "qmc_data.hpp"

namespace cthyb {

  using namespace triqs::gfs;

  // Measure Legendre Green's function (all blocks)
  struct measure_G_l {

    public:
    measure_G_l(std::optional<G_l_t> &G_l_opt, qmc_data const &data, int n_l, gf_struct_t const &gf_struct);
    void accumulate(mc_weight_t s);
    void collect_results(triqs::mpi::communicator const &c);

    private:
    qmc_data const &data;
    mc_weight_t average_sign;
    G_l_t::view_type G_l;
  };

} // namespace cthyb
