/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2026
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

  class move_worm_insert_F {
    qmc_data &data;
    mc_tools::random_generator &rng;
    double worm_eta;
    h_scalar_t new_atomic_weight, new_atomic_reweighting;
    int component_index = -1;
    time_pt tau_Q, tau_cdag;

    public:
    move_worm_insert_F(qmc_data &data, mc_tools::random_generator &rng, double worm_eta);
    mc_weight_t attempt();
    mc_weight_t accept();
    void reject();
  };

  class move_worm_remove_F {
    qmc_data &data;
    mc_tools::random_generator &rng;
    double worm_eta;
    h_scalar_t new_atomic_weight, new_atomic_reweighting;

    public:
    move_worm_remove_F(qmc_data &data, mc_tools::random_generator &rng, double worm_eta);
    mc_weight_t attempt();
    mc_weight_t accept();
    void reject();
  };

  class move_worm_shift_F {
    qmc_data &data;
    mc_tools::random_generator &rng;
    h_scalar_t new_atomic_weight, new_atomic_reweighting;
    bool shift_Q = true;
    time_pt tau_old, tau_new;
    op_desc op_new;

    public:
    move_worm_shift_F(qmc_data &data, mc_tools::random_generator &rng);
    mc_weight_t attempt();
    mc_weight_t accept();
    void reject();
  };

} // namespace triqs_cthyb
