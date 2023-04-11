/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, H. U.R. Strand, N. Wentzell
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

#define STR(x) #x
#define STRINGIZE(x) STR(x)

#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/utility/time_pt.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp> // gf_struct_t
#include <triqs/stat/histograms.hpp>
#include <triqs/atom_diag/atom_diag.hpp>

#include <itertools/itertools.hpp>

#include "config.hpp"

#include <variant>

namespace triqs_cthyb {

  using namespace triqs::gfs;
  using namespace triqs::mesh;
  using namespace triqs::utility;
  using namespace triqs::stat;
  using namespace itertools;

  using atom_diag = triqs::atom_diag::atom_diag<is_h_scalar_complex>;

  using triqs::hilbert_space::gf_struct_t;
  using triqs::utility::time_pt;
  using op_t        = std::pair<time_pt, int>;
  using histo_map_t = std::map<std::string, histogram>;

  using indices_type = triqs::operators::indices_t;
  using gf_struct_t  = triqs::hilbert_space::gf_struct_t;

  // One-particle Green's function types
  using G_tau_t          = block_gf<imtime, matrix_valued>;
  using G_tau_G_target_t = block_gf<imtime, G_target_t>;
  using G_iw_t           = block_gf<imfreq, matrix_valued>;
  using G_l_t            = block_gf<triqs::gfs::legendre, matrix_valued>;

  // Two-particle Green's function types
  using imtime_cube_mesh_t = prod<imtime, imtime, imtime>;
  using G2_tau_t           = block2_gf<imtime_cube_mesh_t, tensor_valued<4>>;

  using imfreq_cube_mesh_t = prod<imfreq, imfreq, imfreq>;
  using G2_iw_t            = block2_gf<imfreq_cube_mesh_t, tensor_valued<4>>;

  using imfreq_legendre_mesh_t = prod<imfreq, triqs::gfs::legendre, triqs::gfs::legendre>;
  using G2_iwll_t              = block2_gf<imfreq_legendre_mesh_t, tensor_valued<4>>;

  enum class G2_channel { PP, PH, AllFermionic }; // G2 sampling channels

  /// Order of block indices for Block2Gf objects
  enum class block_order { AABB, ABBA };

  inline void h5_write(h5::group h5group, std::string name, block_order const &bo) { h5_write(h5group, name, static_cast<int>(bo)); }

  inline void h5_read(h5::group h5group, std::string name, block_order &bo) {
    int idx;
    h5_read(h5group, name, idx);
    bo = static_cast<block_order>(idx);
  }

} // namespace triqs_cthyb

namespace triqs {
  namespace gfs {

    /// Function template for block2_gf initialization
    template <typename Var_t>
    block2_gf<Var_t, tensor_valued<4>> make_block2_gf(Var_t const &m, triqs::hilbert_space::gf_struct_t const &gf_struct,
                                                      triqs_cthyb::block_order order = triqs_cthyb::block_order::AABB) {

      std::vector<std::vector<gf<Var_t, tensor_valued<4>>>> gf_vecvec;
      std::vector<std::string> block_names;

      for (auto const &[bl1, bl1_size] : gf_struct) {
        block_names.push_back(bl1);
        std::vector<std::string> indices1;
        for (auto idx : range(bl1_size)) indices1.push_back(std::to_string(idx));

        std::vector<gf<Var_t, tensor_valued<4>>> gf_vec;
        for (auto const &[bl2, bl2_size] : gf_struct) {
          std::vector<std::string> indices2;
          for (auto idx : range(bl2_size)) indices2.push_back(std::to_string(idx));
          auto I = std::vector<std::vector<std::string>>{indices1, indices1, indices2, indices2};
          switch (order) {
            case triqs_cthyb::block_order::AABB: gf_vec.emplace_back(m, make_shape(bl1_size, bl1_size, bl2_size, bl2_size), I); break;
            case triqs_cthyb::block_order::ABBA: gf_vec.emplace_back(m, make_shape(bl1_size, bl2_size, bl2_size, bl1_size), I); break;
          }
        }
        gf_vecvec.emplace_back(std::move(gf_vec));
      }
      return make_block2_gf(block_names, block_names, std::move(gf_vecvec));
    }

  } // namespace gfs
} // namespace triqs
