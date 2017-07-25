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

#include "../types.hpp"

namespace cthyb {

  // --------------------------------------------------------------------------
  /// Two-particle Green's function block-measure
  /// specifying the two connected blocks
  /// and the resulting target space
  class G2_measure_t {
    public:
    using target_shape_t = mini_vector<int, 4>;
    struct block_t {
      int idx;
      std::string name;
    };
    block_t b1, b2;
    target_shape_t target_shape;
    G2_measure_t(block_t b1, block_t b2, target_shape_t target_shape) : b1(b1), b2(b2), target_shape(target_shape) {}
  };

  // --------------------------------------------------------------------------
  /// Helper class that keeps track of what blocks to sample
  /// passed on to the g4 measurements
  class G2_measures_t {

    private:
    std::vector<G2_measure_t> measures;

    public:
    const gf_struct_t gf_struct;
    const solve_parameters_t params;

    const std::vector<G2_measure_t> &operator()() { return measures; }

    /// the constructor mangles the parameters, especially params.measure_G2_blocks
    /// and populates the std::vector<g4_measure_t> measures
    G2_measures_t(const G_tau_t &_Delta_tau, const gf_struct_t &gf_struct, const solve_parameters_t &params) : gf_struct(gf_struct), params(params) {

      auto G2_blocks_to_measure = params.measure_G2_blocks;

      // Measure all blocks
      if (G2_blocks_to_measure.empty()) {
        for (auto const &bn1 : gf_struct) {
          for (auto const &bn2 : gf_struct) { G2_blocks_to_measure.emplace(bn1.first, bn2.first); }
        }
      } else { // Check the blocks we've been asked to measure
        for (auto const &bn : G2_blocks_to_measure) {
          if (!gf_struct.count(bn.first)) TRIQS_RUNTIME_ERROR << "Invalid left block name " << bn.first << " for G^2 measurement";
          if (!gf_struct.count(bn.second)) TRIQS_RUNTIME_ERROR << "Invalid right block name " << bn.second << " for G^2 measurement";
        }
      }

      std::map<std::string, int> block_name_to_index;
      auto block_names = _Delta_tau.block_names();
      for (auto block_idx : range(block_names.size())) {
        std::string block_name          = block_names[block_idx];
        block_name_to_index[block_name] = block_idx;
      }

      for (auto const &block_pair : G2_blocks_to_measure) {

        auto const &bn1 = std::get<0>(block_pair);
        auto const &bn2 = std::get<1>(block_pair);

        auto b1 = block_name_to_index[bn1];
        auto b2 = block_name_to_index[bn2];

        int s1 = _Delta_tau[b1].target_shape()[0];
        int s3 = _Delta_tau[b2].target_shape()[0];
        int s2 = params.measure_G2_block_order == AABB ? s1 : s3;
        int s4 = params.measure_G2_block_order == AABB ? s3 : s1;

        mini_vector<int, 4> target_shape{s1, s2, s3, s4};

        G2_measure_t measure{{b1, bn1}, {b2, bn2}, target_shape};

        measures.push_back(measure);
      }
    }
  };

  // --------------------------------------------------------------------------
} // namespace cthyb
