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
}
