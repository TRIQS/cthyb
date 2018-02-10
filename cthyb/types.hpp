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

#define STR(x) #x
#define STRINGIZE(x) STR(x)

#include <triqs/gfs.hpp>
#include <triqs/utility/time_pt.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp> // gf_struct_t
#include <triqs/statistics/histograms.hpp>
#include <triqs/atom_diag/atom_diag.hpp>

#include "config.hpp"
#include <variant>

namespace std {
  inline string to_string(string const &str) { return str; }
  inline string to_string(variant<int, string> const &var_int_string) {
    stringstream ss;
    ss << var_int_string;
    return ss.str();
  }
} // namespace std

namespace cthyb {

  using namespace triqs::gfs;
  using namespace triqs::utility;
  using namespace triqs::statistics;

  using atom_diag = triqs::atom_diag::atom_diag<is_h_scalar_complex>;

  using triqs::utility::time_pt;
  using op_t        = std::pair<time_pt, int>;
  using histo_map_t = std::map<std::string, histogram>;

  using indices_type = triqs::operators::indices_t;
  using gf_struct_t  = triqs::hilbert_space::gf_struct_t;

  // One-particle Green's function types
  using G_tau_t          = block_gf<imtime, matrix_valued>;
  using G_tau_G_target_t = block_gf<imtime, G_target_t>;
  using G_iw_t           = block_gf<imfreq, matrix_valued>;
  using G_l_t            = block_gf<legendre, matrix_valued>;

  // Two-particle Green's function types
  using imtime_cube_mesh_t = cartesian_product<imtime, imtime, imtime>;
  using G2_tau_t           = block2_gf<imtime_cube_mesh_t, tensor_valued<4>>;

  using imfreq_cube_mesh_t = cartesian_product<imfreq, imfreq, imfreq>;
  using G2_iw_t            = block2_gf<imfreq_cube_mesh_t, tensor_valued<4>>;

  using imfreq_legendre_mesh_t = cartesian_product<imfreq, legendre, legendre>;
  using G2_iwll_t              = block2_gf<imfreq_legendre_mesh_t, tensor_valued<4>>;

  enum class G2_channel { PP, PH, AllFermionic }; // G2 sampling channels

  /// Order of block indices for Block2Gf objects
  enum class block_order { AABB, ABBA };
  
  inline void h5_write(triqs::h5::group h5group, std::string name, block_order const &bo) {
    h5_write(h5group, name, static_cast<int>(bo));
  }

  inline void h5_read(triqs::h5::group h5group, std::string name, block_order &bo) {
    int idx;
    h5_read(h5group, name, idx);
    bo = static_cast<block_order>(idx);
  }

} // namespace cthyb

namespace triqs {
  namespace gfs {

    /// The structure of the gf : block_name -> [...]= list of indices (int/string). FIXME Change to pair of vec<str> and vec<int> or vec<pair<str,int>>
    using block_gf_structure_t = std::map<std::string, std::vector<std::variant<int, std::string>>>;
    /// Function template for block_gf initialization
    template <typename Val_t, typename Var_t> block_gf<Var_t, Val_t> make_block_gf(gf_mesh<Var_t> const &m, block_gf_structure_t const &gf_struct) {

      std::vector<gf<Var_t, Val_t>> gf_vec;
      std::vector<std::string> block_names;

      //for (auto const & [ bname, idx_lst ] : gf_struct) { // C++17
      for (auto const &bl : gf_struct) {
        auto &bname  = bl.first;
        auto bl_size = bl.second.size();
        block_names.push_back(bname);
        std::vector<std::string> indices;
        for (auto const &var : bl.second) visit([&indices](auto &&arg) { indices.push_back(std::to_string(arg)); }, var);
        gf_vec.emplace_back(m, make_shape(bl_size, bl_size), std::vector<std::vector<std::string>>{indices, indices});
      }

      return make_block_gf(std::move(block_names), std::move(gf_vec));
    }

    // default to matrix_valued gf
    template <typename Var_t> block_gf<Var_t, matrix_valued> make_block_gf(gf_mesh<Var_t> const &m, block_gf_structure_t const &gf_struct) {
      return make_block_gf<matrix_valued, Var_t>(m, gf_struct);
    }

    /// Function template for block2_gf initialization
    template <typename Var_t>
    block2_gf<Var_t, tensor_valued<4>> make_block2_gf(gf_mesh<Var_t> const &m, block_gf_structure_t const &gf_struct,
                                                      cthyb::block_order order = cthyb::block_order::AABB) {

      std::vector<std::vector<gf<Var_t, tensor_valued<4>>>> gf_vecvec;
      std::vector<std::string> block_names;

      for (auto const &bl1 : gf_struct) {
        auto &bname  = bl1.first;
        int bl1_size = bl1.second.size();
        block_names.push_back(bname);
        std::vector<std::string> indices1;
        for (auto const &var : bl1.second) visit([&indices1](auto &&arg) { indices1.push_back(std::to_string(arg)); }, var);

        std::vector<gf<Var_t, tensor_valued<4>>> gf_vec;
        for (auto const &bl2 : gf_struct) {
          int bl2_size = bl2.second.size();
          std::vector<std::string> indices2;
          for (auto const &var : bl2.second) visit([&indices2](auto &&arg) { indices2.push_back(std::to_string(arg)); }, var);
          auto I = std::vector<std::vector<std::string>>{indices1, indices1, indices2, indices2};
          switch (order) {
            case cthyb::block_order::AABB: gf_vec.emplace_back(m, make_shape(bl1_size, bl1_size, bl2_size, bl2_size), I); break;
            case cthyb::block_order::ABBA: gf_vec.emplace_back(m, make_shape(bl1_size, bl2_size, bl2_size, bl1_size), I); break;
          }
        }
        gf_vecvec.emplace_back(std::move(gf_vec));
      }
      return make_block2_gf(block_names, block_names, std::move(gf_vecvec));
    }

  } // namespace gfs
} // namespace triqs
