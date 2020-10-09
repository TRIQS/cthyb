/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, N. Wentzell, I. Krivenko
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

#include <cmath>
#include <array>
#include <numeric>
#include <utility>
#include <triqs/arrays.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/utility/tuple_tools.hpp>

#include <nfft3.h>

namespace triqs {
  namespace experimental {

    using namespace triqs::arrays;
    using namespace triqs::gfs;
    using namespace triqs::mesh;

    using dcomplex = std::complex<double>;

    // NFFT buffer
    template <int Rank> class nfft_buf_t {

      template <typename = std::make_index_sequence<Rank>> struct imfreq_product;
      template <std::size_t... Is> struct imfreq_product<std::index_sequence<Is...>> { using type = prod<decltype(Is, imfreq{})...>; };

      public:
      using freq_mesh_t = gf_mesh<typename imfreq_product<>::type>;

      nfft_buf_t(freq_mesh_t const &fiw_mesh, array_view<dcomplex, Rank> fiw_arr, int buf_size, bool do_checks = false)
         : fiw_mesh(fiw_mesh), fiw_arr(fiw_arr), buf_size(buf_size), do_checks(do_checks), plan_ptr(std::make_unique<nfft_plan>()), buf_counter(0) {

        // Initialize NFFT plans
        all_fermion   = true;
        common_factor = 1;
        triqs::tuple::for_each_enumerate(fiw_mesh.components(), [this](int r, gf_mesh<imfreq> const &m) {
          if (m.domain().statistic == Fermion) {
            index_shifts[r] = 0;
            common_factor *= (m.size() / 2) % 2 ? -1 : 1;
          } else {
            all_fermion     = false;
            index_shifts[r] = 1; // For bosons we discard the most negative frequency
            common_factor *= ((m.size() - 1) / 2) % 2 ? -1 : 1;
            if (m.size() < 5) {
              std::cerr << " ERROR: nfft_buf_t needs more bosonic frequencies.\n";
              exit(0);
            }
          }
        });
        std::array<int, Rank> buf_extents = fiw_mesh.size_of_components() + index_shifts;

        if (!all_fermion) nfft_indexmap = indexmaps::cuboid::domain_t<Rank>(buf_extents);

        // -- Default init
        //nfft_init(plan_ptr.get(), Rank, buf_extents.ptr(), buf_size);

        /// compute the next highest power of 2 of 32-bit v
        auto next_power_of_two = [](unsigned int v) {
          v--;
          v |= v >> 1;
          v |= v >> 2;
          v |= v >> 4;
          v |= v >> 8;
          v |= v >> 16;
          v++;
          return v;
        };

        // Init nfft_plan
        std::array<int, Rank> extents_fftw;
        for (int i = 0; i < Rank; i++) extents_fftw[i] = 2 * next_power_of_two(buf_extents[i]);

        unsigned nfft_flags = PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE | NFFT_SORT_NODES;
        unsigned fftw_flags = FFTW_ESTIMATE | FFTW_DESTROY_INPUT;

        int m = 6; // Truncation order for the window functions
        nfft_init_guru(plan_ptr.get(), Rank, buf_extents.ptr(), buf_size, extents_fftw.ptr(), m, nfft_flags, fftw_flags);
      }

      ~nfft_buf_t() {
        if (buf_counter != 0) std::cerr << " WARNING: Points in NFFT Buffer lost \n";
        if (plan_ptr) nfft_finalize(plan_ptr.get());
      }

      // nfft_buffer needs to be uncopyable, because nfft_plan contains raw pointers
      nfft_buf_t(nfft_buf_t const &) = delete;
      nfft_buf_t(nfft_buf_t &&)      = default;
      nfft_buf_t &operator=(nfft_buf_t const &) = delete;
      nfft_buf_t &operator=(nfft_buf_t &&) = default;

      /// Rebind nfft buffer to new accumulation container of same shape
      void rebind(array_view<dcomplex, Rank> new_fiw_arr) {
        flush();
        assert(get_shape(new_fiw_arr) == get_shape(fiw_arr) || get_shape(fiw_arr) == get_shape(array_view<dcomplex, Rank>{}));
        fiw_arr.rebind(new_fiw_arr);
      }

      // Insert tau-vector {\tau_1, \tau_2, ... } \in [0,\beta_1)x[0,\beta_2)x...
      // and the corresponding f(tau) into the NFFT buffer
      void push_back(std::array<double, Rank> const &tau_arr, dcomplex ftau) {

        assert(plan_ptr);

        // Write the set of shifted and normalized tau values (i. e. x values) to
        // the NFFT buffer and sum \tau/\beta for fermions
        double sum_tau_beta = 0;
        triqs::tuple::for_each_enumerate(fiw_mesh.components(), [&tau_arr, &sum_tau_beta, this](int r, gf_mesh<imfreq> const &m) {
          double tau  = tau_arr[r];
          double beta = m.domain().beta;

          // Note: Nfft multi-arrays are stored in flattened arrays (c-order)
          x_arr()[buf_counter * Rank + r] = tau_arr[r] / beta - 0.5; // \in [-0.5, 0.5)

          if (m.domain().statistic == Fermion) sum_tau_beta += tau / beta;
        });

        // Write f(x) to nfft_plan-> The prefactor accounts for the Pi/beta offset
        // in fermionic Matsubaras
        fx_arr()[buf_counter] = std::exp(1i * M_PI * sum_tau_beta) * ftau;

        ++buf_counter;

        // If buffer is full, perform transform
        if (is_full()) {
          do_nfft();
          buf_counter = 0;
        }
      }

      /// Flush contents of the nfft buffer
      void flush() {

        assert(plan_ptr);

        if (is_empty()) return;

        // Trivial initialization of the remaining points
        for (int i = buf_counter; i < buf_size; ++i) {
          fx_arr()[i] = 0.0;
          for (int r = 0; r < Rank; ++r) x_arr()[i * Rank + r] = -0.5 + double(i) / buf_size;
        }
        do_nfft();
        buf_counter = 0;
      }

      // Function to check whether buffer is empty
      bool is_empty() const { return buf_counter == 0; }

      // Function to check whether buffer is filled
      bool is_full() const { return buf_counter >= buf_size; }

      private:
      // Imaginary frequency (multi-)mesh
      freq_mesh_t fiw_mesh;

      // TRIQS array to contain the final NFFT output in matsubara frequencies
      array_view<dcomplex, Rank> fiw_arr;

      // Are all mesh components fermionic?
      bool all_fermion;

      // Index map for plan_ptr->f_hat
      indexmaps::cuboid::map<Rank> nfft_indexmap;

      // Bosonic components of fiw_mesh indices must be shifted by 1,
      // since we want to ignore the most negative frequencies in plan_ptr->f_hat
      std::array<std::size_t, Rank> index_shifts;

      // Common prefactor for the transformation result
      int common_factor;

      // Number of tau points for the nfft
      int buf_size;

      // Switch for testing in nfft
      bool do_checks;

      // Pointer to NFFT plan
      std::unique_ptr<nfft_plan> plan_ptr;

      // Counter for elements currently in the buffer
      int buf_counter;

      // Get pointer to array containing x values for the NFFT transform
      double *x_arr() { return plan_ptr->x; }

      // Get pointer to array containing f(x) values for the NFFT transform
      dcomplex *fx_arr() { return reinterpret_cast<dcomplex *>(plan_ptr->f); }

      // Get pointer to array containing the NFFT output h(k)
      const dcomplex *fk_arr() const { return reinterpret_cast<dcomplex *>(plan_ptr->f_hat); }

      // Perform NFFT transform and accumulate inside fiw_arr
      void do_nfft() {

        // nfft_adjoint() uses a window function (Kaiser-Bessel by default)
        // that cannot be constructed for plan_ptr->N[i] < plan_ptr->m.
        // In the small N[i] case one has to call nfft_adjoint_direct() instead,
        // which is also faster for the smaller N[i].
        //
        // C.f. https://github.com/NFFT/nfft/issues/34
        auto N_min = *std::min_element(plan_ptr->N, plan_ptr->N + plan_ptr->d);
        if (N_min <= plan_ptr->m) {

          // Execute transform
          nfft_adjoint_direct(plan_ptr.get());

        } else {

// NFFT Library precomputation and checks
#ifdef NFFT_OLD_API
          if (plan_ptr->nfft_flags & PRE_ONE_PSI) nfft_precompute_one_psi(plan_ptr.get());
#else
          if (plan_ptr->flags & PRE_ONE_PSI) nfft_precompute_one_psi(plan_ptr.get());
#endif
          if (do_checks) { // Check validity of NFFT parameters
            const char *error_str = nfft_check(plan_ptr.get());
            if (error_str != 0) TRIQS_RUNTIME_ERROR << "Error in NFFT module: " << error_str << "\n";
          }

          // Execute transform
          nfft_adjoint(plan_ptr.get());
        }

        // Accumulate results in fiw_arr. Care to normalize results afterwards
        if (all_fermion) {
          int count = 0;
          for (auto it = fiw_arr.begin(); it != fiw_arr.end(); ++it, ++count) {
            int factor = (sum(it.indices()) % 2 ? -1 : 1);
            *it += fk_arr()[count] * factor * common_factor;
          }
        } else {
          for (auto it = fiw_arr.begin(); it != fiw_arr.end(); ++it) {
            int count  = nfft_indexmap[it.indices() + index_shifts];
            int factor = (sum(it.indices()) % 2 ? -1 : 1);
            *it += fk_arr()[count] * factor * common_factor;
          }
        }
      }
    };
  }
}
