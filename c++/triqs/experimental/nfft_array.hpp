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
    using triqs::utility::time_pt;

    // NFFT transform of an array-valued function of MeshRank tau arguments
    template <int MeshRank, int TargetRank> class nfft_array_t {

      public:
      using freq_mesh_t = typename nfft_buf_t<MeshRank>::freq_mesh_t;
      using res_gf_t    = gf<typename freq_mesh_t::var_t, tensor_valued<TargetRank>>;

      nfft_array_t() = default;

      // fiw_mesh - Matsubara frequency mesh
      // shape - target shape
      nfft_array_t(freq_mesh_t const &fiw_mesh, mini_vector<int, TargetRank> const &shape) : indexmap(shape), max_n_tau(1), result(fiw_mesh, shape) {
        long n = result.data().domain().number_of_elements();
        buffers.reserve(n);
        for (int i : range(n))
#ifdef NDEBUG
          buffers.emplace_back(result.mesh(), max_n_tau, false);
#else
          buffers.emplace_back(result.mesh(), max_n_tau, true);
#endif
      }

      // Resize all NFFT buffers if their capacity is insufficient
      void resize_bufs(int n_tau) {
        if (n_tau > max_n_tau) {
          max_n_tau = n_tau;
          for (auto &buf : buffers)
#ifdef NDEBUG
            buf = {result.mesh(), max_n_tau, false};
#else
            buf = {result.mesh(), max_n_tau, true};
#endif
        }
      }

      // Add a new element to the NFFT buffer
      void push_back(std::array<double, MeshRank> const &tau_arr, mini_vector<int, TargetRank> const &ind_arr, dcomplex fxy) {
        select_buffer(ind_arr, std14::make_index_sequence<TargetRank>()).push_back(tau_arr, fxy);
      }

      // Run transformation
      void transform() {
        for (auto const &ind_arr : indexmap.domain()) { transform_impl(ind_arr, std14::make_index_sequence<TargetRank>()); }
      }

      // Access the result g_{ab...}(iw_1, iw_2, ...)
      res_gf_t const &operator()() { return result; }

      private:
      template <size_t... Is> inline nfft_buf_t<MeshRank> &select_buffer(mini_vector<int, TargetRank> const &ind_arr, std14::index_sequence<Is...>) {
        return buffers[indexmap(ind_arr[Is]...)];
      }

      template <size_t... Is> inline void transform_impl(mini_vector<int, TargetRank> const &ind_arr, std14::index_sequence<Is...>) {
        auto &buf = select_buffer(ind_arr, std14::make_index_sequence<TargetRank>());
        if (buf.is_empty()) {
          result.data()(ellipsis(), ind_arr[Is]...) = 0;
        } else {
          buf.flush();
          buf.fill_array(result.data()(ellipsis(), ind_arr[Is]...));
        }
      }

      // Index map for the target array
      indexmaps::cuboid::map<TargetRank> indexmap;

      // Maximum number of tau-pairs over all elements
      int max_n_tau;

      // NFFT buffers
      std::vector<nfft_buf_t<MeshRank>> buffers;

      // NFFT transformation result
      res_gf_t result;
    };
  }
}
