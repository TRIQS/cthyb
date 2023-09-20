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
#include <triqs/utility/time_pt.hpp>
#include <triqs/experimental/nfft_buf.hpp>

namespace triqs {
  namespace experimental {

    using namespace nda;
    using namespace triqs::gfs;
    using namespace triqs::mesh;
    using triqs::utility::time_pt;

    // NFFT transform of an array-valued function of MeshRank tau arguments
    template <int MeshRank, int TargetRank> class nfft_array_t {

      public:
      using freq_mesh_t                = typename nfft_buf_t<MeshRank>::freq_mesh_t;
      using res_gf_t                   = gf<freq_mesh_t, tensor_valued<TargetRank>>;
      static constexpr int result_rank = MeshRank + TargetRank;

      nfft_array_t() = default;

      nfft_array_t(nfft_array_t const &other) = delete;
      nfft_array_t(nfft_array_t &&other) : indexmap(other.indexmap), buffers(std::move(other.buffers)) { fiw_arr.rebind(other.fiw_arr); }

      nfft_array_t &operator=(nfft_array_t const &other) = delete;
      nfft_array_t &operator=(nfft_array_t &&other) {
        indexmap = other.indexmap;
        buffers  = std::move(other.buffers);
        fiw_arr.rebind(other.fiw_arr);
        return *this;
      }

      // fiw_mesh - Matsubara frequency mesh
      // fiw_arr_ - array to contain the final NFFT output
      // buf_sizes - sizes of NFFT buffers
      nfft_array_t(freq_mesh_t const &fiw_mesh, array_view<dcomplex, result_rank> fiw_arr_, array<long, TargetRank> const &buf_sizes)
         : indexmap(make_target_shape(fiw_arr_.shape())), fiw_arr(fiw_arr_) {
        buffers.reserve(indexmap.size());
        for_each(buf_sizes.shape(), [this, &fiw_mesh, &buf_sizes](auto... ind) {
#ifdef NDEBUG
          buffers.emplace_back(fiw_mesh, fiw_arr(ellipsis(), ind...), buf_sizes(ind...), false);
#else
          buffers.emplace_back(fiw_mesh, fiw_arr(ellipsis(), ind...), buf_sizes(ind...), true);
#endif
        });
      }

      // Add a new element to the NFFT buffer
      void push_back(std::array<double, MeshRank> const &tau_arr, std::array<int, TargetRank> const &ind_arr, dcomplex fxy) {
        select_buffer(ind_arr, std::make_index_sequence<TargetRank>()).push_back(tau_arr, fxy);
      }

      // Run transformation
      void flush() {
        for (auto &buf : buffers) buf.flush();
      }

      private:
      template <size_t... Is> inline nfft_buf_t<MeshRank> &select_buffer(std::array<int, TargetRank> const &ind_arr, std::index_sequence<Is...>) {
        return buffers[indexmap(ind_arr[Is]...)];
      }

      std::array<long, TargetRank> make_target_shape(std::array<long, result_rank> const &shape) {
        std::array<long, TargetRank> res;
        for (int n = 0; n < TargetRank; ++n) res[n] = shape[n + MeshRank];
        return std::array<long, TargetRank>(res);
      }

      // Index map for the target array
      nda::idx_map<TargetRank, 0, C_stride_order<TargetRank>, layout_prop_e::none> indexmap;

      // NFFT buffers
      std::vector<nfft_buf_t<MeshRank>> buffers;

      // NFFT transformation result
      array_view<dcomplex, result_rank> fiw_arr;
    };
  } // namespace experimental
} // namespace triqs
