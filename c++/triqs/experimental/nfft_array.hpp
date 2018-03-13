/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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

#include <vector>
#include <triqs/clef.hpp>
#include <triqs/utility/time_pt.hpp>
#include <triqs/experimental/nfft_buf.hpp>

namespace triqs {
  namespace experimental {

    using namespace triqs::arrays;
    using namespace triqs::gfs;
    using triqs::utility::time_pt;

    // NFFT transform of an array-valued function of MeshRank tau arguments
    template <int MeshRank, int TargetRank> class nfft_array_t {

      public:
      using freq_mesh_t                = typename nfft_buf_t<MeshRank>::freq_mesh_t;
      using res_gf_t                   = gf<typename freq_mesh_t::var_t, tensor_valued<TargetRank>>;
      static constexpr int result_rank = MeshRank + TargetRank;

      nfft_array_t() = default;

      // fiw_mesh - Matsubara frequency mesh
      // fiw_arr_ - array to contain the final NFFT output
      // buf_sizes - sizes of NFFT buffers
      nfft_array_t(freq_mesh_t const &fiw_mesh, array_view<dcomplex, result_rank> fiw_arr_, array<int, TargetRank> const &buf_sizes)
         : indexmap(make_target_shape(fiw_arr_.shape())), fiw_arr(fiw_arr_) {
        buffers.reserve(indexmap.domain().number_of_elements());
        foreach (buf_sizes, [this, &fiw_mesh, &buf_sizes](auto... ind) {
#ifdef NDEBUG
          buffers.emplace_back(fiw_mesh, fiw_arr(ellipsis(), ind...), buf_sizes(ind...), false);
#else
          buffers.emplace_back(fiw_mesh, fiw_arr(ellipsis(), ind...), buf_sizes(ind...), true);
#endif
        })
          ;
      }

      // Add a new element to the NFFT buffer
      void push_back(std::array<double, MeshRank> const &tau_arr, mini_vector<int, TargetRank> const &ind_arr, dcomplex fxy) {
        select_buffer(ind_arr, std14::make_index_sequence<TargetRank>()).push_back(tau_arr, fxy);
      }

      // Run transformation
      void flush() {
        for (auto &buf : buffers) buf.flush();
      }

      private:
      template <size_t... Is> inline nfft_buf_t<MeshRank> &select_buffer(mini_vector<int, TargetRank> const &ind_arr, std14::index_sequence<Is...>) {
        return buffers[indexmap(ind_arr[Is]...)];
      }

      mini_vector<int, TargetRank> make_target_shape(mini_vector<int, result_rank> const &shape) {
        std::vector<int> res(TargetRank);
        for (int n = 0; n < TargetRank; ++n) res[n] = shape[n + MeshRank];
        return mini_vector<int, TargetRank>(res);
      }

      // Index map for the target array
      indexmaps::cuboid::map<TargetRank> indexmap;

      // NFFT buffers
      std::vector<nfft_buf_t<MeshRank>> buffers;

      // NFFT transformation result
      array_view<dcomplex, result_rank> fiw_arr;
    };
  } // namespace experimental
} // namespace triqs
