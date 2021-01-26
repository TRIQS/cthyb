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

  // Move a C or C^dagger operator to a different time
  class move_shift_operator {

    qmc_data &data;
    configuration &config;
    mc_tools::random_generator &rng;
    histogram *histo_proposed, *histo_accepted; // Analysis histograms
    double dtau;
    h_scalar_t new_atomic_weight, new_atomic_reweighting;
    time_pt tau_old, tau_new;
    op_desc op_old, op_new;
    using det_type = det_manip::det_manip_basic<qmc_data::delta_block_adaptor>;
    det_type::RollDirection roll_direction;
    int block_index;

    histogram *add_histo(std::string const &name, histo_map_t *histos);

    public:
    move_shift_operator(qmc_data &data, mc_tools::random_generator &rng, histo_map_t *histos);
    mc_weight_t attempt();
    mc_weight_t accept();
    void reject();
  };
}
