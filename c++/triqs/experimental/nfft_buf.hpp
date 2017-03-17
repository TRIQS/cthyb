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

#include <array>
#include <numeric>
#include <utility>
#include <triqs/arrays.hpp>
#include <triqs/gfs.hpp>
#include <triqs/utility/tuple_tools.hpp>

#define NFFT_PRECISION_DOUBLE
#ifdef DNFFT_OLD_API
#include <nfft3.h>
#else
#include <nfft3mp.h>
#endif

namespace triqs {
  namespace experimental {

    // FIXME: remove when merged into TRIQS library master
    template <typename T, int Rank> T sum(triqs::arrays::mini_vector<T, Rank> const &a) {
      T res = {};
      for (int i = 0; i < Rank; ++i) res += a[i];
      return res;
    }

    using namespace triqs::arrays;
    using namespace triqs::gfs;

    using dcomplex = std::complex<double>;

    // NFFT buffer
    template <int Rank> class nfft_buf_t {

      template <typename = std::make_index_sequence<Rank>> struct imfreq_product;
      template <std::size_t... Is> struct imfreq_product<std::index_sequence<Is...>> { using type = cartesian_product<decltype(Is, imfreq{})...>; };

      public:
      using freq_mesh_t = gf_mesh<typename imfreq_product<>::type>;

      nfft_buf_t(freq_mesh_t const &fiw_mesh, int buf_size, bool do_checks = false)
         : fiw_mesh(fiw_mesh), plan_ptr(std::make_unique<nfft_plan>()), buf_counter(0), buf_size(buf_size), do_checks(do_checks) {

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
          }
        });
        mini_vector<int, Rank> buf_extents = fiw_mesh.size_of_components() + index_shifts;

        if (!all_fermion) nfft_indexmap = indexmaps::cuboid::domain_t<Rank>(buf_extents);

        nfft_init(plan_ptr.get(), Rank, buf_extents.ptr(), buf_size);
      }

      ~nfft_buf_t() {
        if (plan_ptr) nfft_finalize(plan_ptr.get());
      }

      // nfft_buffer needs to be uncopyable, because nfft_plan contains raw pointers
      nfft_buf_t(nfft_buf_t const &) = delete;
      nfft_buf_t(nfft_buf_t &&)      = default;
      nfft_buf_t &operator=(nfft_buf_t const &) = delete;
      nfft_buf_t &operator=(nfft_buf_t &&) = default;

      // Insert tau-vector {\tau_1, \tau_2, ... } \in [0,\beta_1)x[0,\beta_2)x...
      // and the corresponding f(tau) into the NFFT buffer
      void push_back(std::array<double, Rank> const &tau_arr, dcomplex ftau) {

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
        fx_arr()[buf_counter] = std::exp(1_j * M_PI * sum_tau_beta) * ftau;

        ++buf_counter;

        // If buffer is full, perform transform
        if (is_full()) {
          do_nfft();
          buf_counter = 0;
        }
      }

      /// Flush contents of the nfft buffer
      void flush() {
        // Trivial initialization of the remaining points
        for (int i = buf_counter; i < buf_size; ++i) {
          fx_arr()[i] = 0.0;
          for (int r = 0; r < Rank; ++r) x_arr()[i * Rank + r] = -0.5 + double(i) / buf_size;
        }
        do_nfft();
        buf_counter = 0;
      }

      void fill_array(array_view<dcomplex, Rank> data) {
        if (all_fermion) {
          int count = 0;
          for (auto it = data.begin(); it != data.end(); ++it, ++count) {
            int factor = (sum(it.indices()) % 2 ? -1 : 1);
            *it        = fk_arr()[count] * factor * common_factor;
          }
        } else {
          for (auto it = data.begin(); it != data.end(); ++it) {
            long count = nfft_indexmap[it.indices() + index_shifts];
            int factor = (sum(it.indices()) % 2 ? -1 : 1);
            *it        = fk_arr()[count] * factor * common_factor;
          }
        }
      }

      // Function to check whether buffer is empty
      bool is_empty() const { return buf_counter == 0; }

      // Function to check whether buffer is filled
      bool is_full() const { return buf_counter >= buf_size; }

      private:
      // Imaginary frequency (multi-)mesh
      freq_mesh_t fiw_mesh;

      // Are all mesh components fermionic?
      bool all_fermion;

      // Index map for plan_ptr->f_hat
      indexmaps::cuboid::map<Rank> nfft_indexmap;

      // Bosonic components of fiw_mesh indices must be shifted by 1,
      // since we want to ignore the most negative frequencies in plan_ptr->f_hat
      mini_vector<std::size_t, Rank> index_shifts;

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
    };
  }
}
