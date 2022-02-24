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
#include <triqs/mc_tools.hpp>
#include "../qmc_data.hpp"

namespace triqs_cthyb {

  // Insertion of C, C^dagger operator
  class move_insert_c_cdag {

    qmc_data &data;
    configuration &config;
    mc_tools::random_generator &rng;
    int block_index, block_size;
    histogram *histo_proposed, *histo_accepted; // Analysis histograms
    double dtau;
    h_scalar_t new_atomic_weight, new_atomic_reweighting;
    time_pt tau1, tau2;
    op_desc op1, op2;

    public:
    move_insert_c_cdag(int block_index, int block_size, std::string const &block_name, qmc_data &data, mc_tools::random_generator &rng,
                       histo_map_t *histos);

    mc_weight_t attempt();
    mc_weight_t accept();
    void reject();
  };
}
