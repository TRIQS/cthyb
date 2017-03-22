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

namespace cthyb {

    using namespace triqs::arrays;
    using namespace triqs::gfs;
    using triqs::utility::time_pt;
    using triqs::experimental::nfft_buf_t;

    // NFFT transform of an array-valued function of MeshRank tau arguments
    template <int MeshRank, int TargetRank> class nfft_array_t {

      public:
      using freq_mesh_t = typename nfft_buf_t<MeshRank>::freq_mesh_t;
      using res_gf_t    = gf<typename freq_mesh_t::var_t, tensor_valued<TargetRank>>;

      nfft_array_t() = default;

      // fiw_mesh - Matsubara frequency mesh
      // shape - target shape
      nfft_array_t(freq_mesh_t const &fiw_mesh, mini_vector<int, TargetRank> const &shape, array<int, TargetRank> const &buf_sizes) :
      // FORTRAN_LAYOUT ensures that every nfft_buf receives a contiguous view
      indexmap(shape), result(fiw_mesh, shape, FORTRAN_LAYOUT) {
        auto & data = result.data();
        long n = data.domain().number_of_elements();
        buffers.reserve(n);
        foreach(buf_sizes, [this, &data, &buf_sizes](auto... ind) {
#ifdef NDEBUG
          buffers.emplace_back(result.mesh(), data(ellipsis(), ind...), buf_sizes(ind...), false);
#else
          buffers.emplace_back(result.mesh(), data(ellipsis(), ind...), buf_sizes(ind...), true);
#endif
        });
      }

      // Add a new element to the NFFT buffer
      void push_back(std::array<double, MeshRank> const &tau_arr, mini_vector<int, TargetRank> const &ind_arr, dcomplex fxy) {
        select_buffer(ind_arr, std14::make_index_sequence<TargetRank>()).push_back(tau_arr, fxy);
      }

      // Run transformation
      void flush() {
        for (auto & buf : buffers) buf.flush();
      }

      // Access the result g_{ab...}(iw_1, iw_2, ...)
      res_gf_t &operator()() { return result; }

      private:
      template <size_t... Is> inline nfft_buf_t<MeshRank> &select_buffer(mini_vector<int, TargetRank> const &ind_arr, std14::index_sequence<Is...>) {
        return buffers[indexmap(ind_arr[Is]...)];
      }

      // Index map for the target array
      indexmaps::cuboid::map<TargetRank> indexmap;

      // NFFT buffers
      std::vector<nfft_buf_t<MeshRank>> buffers;

      // NFFT transformation result
      res_gf_t result;
    };
}
