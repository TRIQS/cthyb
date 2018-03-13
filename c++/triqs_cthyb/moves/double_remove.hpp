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
#include <algorithm>
#include "../qmc_data.hpp"

namespace triqs_cthyb {

  // Removal of 2 C and 2 C^dagger operator
  class move_remove_c_c_cdag_cdag {

    qmc_data &data;
    configuration &config;
    mc_tools::random_generator &rng;
    int block_index1, block_index2, block_size1, block_size2;
    // Analysis histograms
    histogram *histo_proposed1, *histo_proposed2;
    histogram *histo_accepted1, *histo_accepted2;
    double dtau1, dtau2;
    h_scalar_t new_atomic_weight, new_atomic_reweighting;
    time_pt tau1, tau2, tau3, tau4;

    histogram *add_histo(std::string const &name, histo_map_t *histos);

    public:
    move_remove_c_c_cdag_cdag(int block_index1, int block_index2, int block_size1, int block_size2, std::string const &block_name1,
                              std::string const &block_name2, qmc_data &data, mc_tools::random_generator &rng, histo_map_t *histos);

    mc_weight_t attempt();
    mc_weight_t accept();
    void reject();
  };
}
