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

#include <triqs/experimental/nfft_matrix.hpp>

#ifndef NDEBUG
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#endif

namespace triqs {

  namespace experimental {

#ifdef NDEBUG
    bool nfft_matrix_t::do_checks = false;
#else
    bool nfft_matrix_t::do_checks = true;
#endif

    nfft_matrix_t::nfft_matrix_t(int size1, int size2, double beta, int n_freq1, int n_freq2)
       : size1(size1), size2(size2), max_n_tau(1), result({{{beta, Fermion, n_freq1}, {beta, Fermion, n_freq2}}, {size1, size2}}) {
      buffers.reserve(size1 * size2);
      for (int a : range(size1))
        for (int b : range(size2)) buffers.emplace_back(result.mesh(), max_n_tau, do_checks);
    }

    void nfft_matrix_t::resize_bufs(int n_tau) {
      if (n_tau > max_n_tau) {
        max_n_tau                     = n_tau;
        for (auto &buf : buffers) buf = {result.mesh(), max_n_tau, do_checks};
      }
    }

    void nfft_matrix_t::push_back(std::pair<time_pt, int> const &x, std::pair<time_pt, int> const &y, dcomplex fxy) {
      buffers[x.second * size2 + y.second].push_back({double(x.first), double(y.first)}, fxy);
    }

    auto nfft_matrix_t::operator()() const -> res_gf_t const & { return result; }

    void nfft_matrix_t::transform() {
      for (int a : range(size1))
        for (int b : range(size2)) {
          auto &buf = buffers[a * size2 + b];
          if (buf.is_empty()) {
            slice_target_to_scalar(result, a, b) = 0;
          } else {
            buf.flush();
            buf.fill_gf(slice_target_to_scalar(result, a, b));
          }
        }
    }
  }
}
