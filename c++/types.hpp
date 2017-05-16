/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, H. U.R. Strand
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
#include <triqs/utility/time_pt.hpp>

#include "config.hpp"

namespace cthyb {

  using namespace triqs::gfs;
  using namespace triqs::utility;
  using namespace triqs::statistics;

  using triqs::utility::time_pt;
  using op_t = std::pair<time_pt, int>;
  using histo_map_t  = std::map<std::string, histogram>;
  using indices_type = triqs::operators::indices_t;
  using gf_struct_t = std::map<std::string, indices_type>;
  
  // One-particle Green's function types
  using g_tau_t = block_gf<imtime>;
  using g_tau_opt_t = std::optional<block_gf<imtime>>;
  using g_tau_view_t = block_gf_view<imtime>;

  using g_tau_g_target_t = block_gf<imtime, g_target_t>;
  
  using g_iw_t = block_gf<imfreq>;
  using g_iw_opt_t = std::optional<block_gf<imfreq>>;
  using g_iw_view_t = block_gf_view<imfreq>;
  
  using g_l_t = block_gf<legendre>;
  using g_l_opt_t = std::optional<block_gf<legendre>>;
  using g_l_view_t = block_gf_view<legendre>;

  // Two-particle Green's function types
  using imtime_cube_mesh_t = cartesian_product<imtime, imtime, imtime>;
  using g4_tau_t = block2_gf<imtime_cube_mesh_t, tensor_valued<4>>;
  using g4_tau_view_t = block2_gf_view<imtime_cube_mesh_t, tensor_valued<4>>;

  using imfreq_cube_mesh_t = cartesian_product<imfreq, imfreq, imfreq>;
  using g4_iw_t = block2_gf<imfreq_cube_mesh_t, tensor_valued<4>>;
  using g4_iw_view_t = block2_gf_view<imfreq_cube_mesh_t, tensor_valued<4>>;

  using imfreq_legendre_mesh_t = cartesian_product<imfreq, legendre, legendre>;
  using g4_wll_t = block2_gf<imfreq_legendre_mesh_t, tensor_valued<4>>;
  using g4_wll_view_t = block2_gf_view<imfreq_legendre_mesh_t, tensor_valued<4>>;
  
} // namespace cthyb
