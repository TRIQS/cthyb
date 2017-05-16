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

  enum g2_channel { PP, PH, AllFermionic };

  namespace measures {

    using g2_channel = cthyb::g2_channel;
    
    class make_measure_name_t {
      public:
      make_measure_name_t(const std::string &bn1, const std::string &bn2) : bn1(bn1), bn2(bn2) {}

      std::string bn1, bn2;
      auto block_quadruple(block_order bo) const {
        return bo == AABB ? (" (" + bn1 + "," + bn1 + "," + bn2 + "," + bn2 + ")") : (" (" + bn1 + "," + bn2 + "," + bn2 + "," + bn1 + ")");
      }
      auto channel_repr(g2_channel channel) {
        if (channel == PP)
          return std::string("pp");
        else if (channel == PH)
          return std::string("ph");
        else if (channel == AllFermionic)
          return std::string("AllFermionic");
        else
          return std::string("unknown");
      }
      auto operator()(imtime, block_order bo) { return std::string("G^2 measure, ImTime") + block_quadruple(bo); }
      auto operator()(imfreq, g2_channel channel, block_order bo) {
        return std::string("G^2 measure, Matsubara, ") + channel_repr(channel) + block_quadruple(bo);
      }
      auto operator()(legendre, g2_channel channel, block_order bo) {
        return std::string("G^2 measure, Legendre, ") + channel_repr(channel) + block_quadruple(bo);
      }
    };

    /// Two-particle Green's function block-measure
    /// specifying the two connected blocks
    /// and the resulting target space
    class g4_measure_t {
      public:
      using target_shape_t = mini_vector<int, 4>;
      struct block_t {
        int idx;
        std::string name;
      };
      block_t b1, b2;
      target_shape_t target_shape;
      g4_measure_t(block_t b1, block_t b2, target_shape_t target_shape) : b1(b1), b2(b2), target_shape(target_shape) {}
    };

    class g4_measures_t {
      
      std::vector<g4_measure_t> measures;

      public:
      const std::vector<g4_measure_t> &operator()() { return measures; }
      
      g4_measures_t(const g_tau_t &_Delta_tau, const gf_struct_t &gf_struct, const solve_parameters_t &params) {

        auto g2_blocks_to_measure = params.measure_g2_blocks;

        // Measure all blocks
        if (g2_blocks_to_measure.empty()) {
          for (auto const &bn1 : gf_struct) {
            for (auto const &bn2 : gf_struct) { g2_blocks_to_measure.emplace(bn1.first, bn2.first); }
          }
        } else { // Check the blocks we've been asked to measure
          for (auto const &bn : g2_blocks_to_measure) {
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

        for (auto const &block_pair : g2_blocks_to_measure) {

          auto const &bn1 = std::get<0>(block_pair);
          auto const &bn2 = std::get<1>(block_pair);

          auto b1 = block_name_to_index[bn1];
          auto b2 = block_name_to_index[bn2];

          int s1 = _Delta_tau[b1].target_shape()[0];
          int s3 = _Delta_tau[b2].target_shape()[0];
          int s2 = params.measure_g2_block_order == AABB ? s1 : s3;
          int s4 = params.measure_g2_block_order == AABB ? s3 : s1;

          mini_vector<int, 4> target_shape{s1, s2, s3, s4};

          g4_measure_t measure{{b1, bn1}, {b2, bn2}, target_shape};

          measures.push_back(measure);
        }
      }
    };

  } // namespace measures

} // namespace cthyb
