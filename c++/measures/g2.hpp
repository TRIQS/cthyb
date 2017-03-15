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
#include <triqs/mpi/base.hpp>
#include <triqs/clef.hpp>
#include "../qmc_data.hpp"
#include "../nfft_buf.hpp"

namespace cthyb {

  using namespace triqs::arrays;

  enum g2_channel { PP, PH };

  namespace details {

    // NFFT transform of a matrix-valued function of two tau arguments
    class nfft_matrix_t {

      // Perform NFFT checks
      static bool do_checks;

      public:
      using res_gf_t = gf<cartesian_product<imfreq, imfreq>, matrix_valued>;

      nfft_matrix_t() = default;

      // size1, size2 - matrix dimensions
      // beta - inverse temperature
      // n_freq1, n_freq2 - number of positive Matsubara frequencies
      // to calculate the Fourier transformation on
      nfft_matrix_t(int size1, int size2, double beta, int n_freq1, int n_freq2);

      // Resize all NFFT buffers if their capacity is insufficient
      void resize_bufs(int n_tau);

      // Add a new matrix element to the NFFT buffer
      void push_back(std::pair<time_pt, int> const &x, std::pair<time_pt, int> const &y, dcomplex fxy);

      // Run transformation
      void transform();

      // Access the result g_{ab}(iw_1, i_w2)
      res_gf_t const &operator()() const;

      private:
      // Matrix sizes
      int size1, size2;

      // Maximum number of tau-pairs over all matrix elements int max_n_tau;
      int max_n_tau;

      // NFFT transformation result
      res_gf_t result;

      // NFFT buffers
      std::vector<nfft_buf_t<2>> buffers;
    };
  }

  // Measure one block of G^2(iomega,inu,inu')
  template <g2_channel Channel, block_order Order> struct measure_g2_inu {

    using g2_view_type = gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>>;
    
    qmc_data const &data;
    g2_view_type g2;
    const int A;           // Index of block A within gf_struct
    const int B;           // Index of block B within gf_struct
    const bool diag_block; // A == B
    mc_weight_t z;
    int64_t num;

    // Objects that perform 2D NFFT transforms
    details::nfft_matrix_t nfft_matrix_ab, nfft_matrix_cd;
    details::nfft_matrix_t nfft_matrix_ad, nfft_matrix_cb;

    measure_g2_inu(int A, int B, g2_view_type g2, qmc_data const &data);
    void accumulate(mc_weight_t s);
    void collect_results(triqs::mpi::communicator const &c);
  };

  // Measure one block of G^2(iomega,l,l')
  template <g2_channel Channel, block_order Order> struct measure_g2_legendre {

    using g2_view_type = gf_view<cartesian_product<imfreq, legendre, legendre>, tensor_valued<4>>;

    qmc_data const &data;
    g2_view_type g2;
    const int A;                 // Index of block A within gf_struct
    const int B;                 // Index of block B within gf_struct
    const bool diag_block;       // A == B
    const size_t size_A, size_B; // Sizes of blocks A and B
    mc_weight_t z;
    int64_t num;

    measure_g2_legendre(int b1, int b2, g2_view_type g2, qmc_data const &data);
    void accumulate(mc_weight_t s);
    void collect_results(triqs::mpi::communicator const &c);
  };
}
