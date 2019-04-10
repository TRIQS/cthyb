/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016, P. Seth, I. Krivenko, H. U.R. Strand, M. Ferrero and O. Parcollet
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

#include "./G2_iw_nfft.hpp"

namespace triqs_cthyb {

  using namespace G2_iw;

  template <G2_channel Channel>
  measure_G2_iw_nfft<Channel>::measure_G2_iw_nfft(std::optional<G2_iw_t> &G2_iw_opt,
                                                  qmc_data const &data,
                                                  G2_measures_t const &G2_measures)
     : measure_G2_iw_base<Channel>(G2_iw_opt, data, G2_measures) {

    // Initialize the nfft_buffers mirroring the matrix M
    {
      M_nfft.resize(M.size());

      std::cout << "G2_iw_nfft NFFT buffer sizes:\n";
      for (auto bidx : range(M.size())) {
        std::string bname = M.block_names()[bidx];

        int buf_size = 10; // default
        if (G2_measures.params.nfft_buf_sizes.count(bname)) {
          buf_size = G2_measures.params.nfft_buf_sizes.at(bname);
        }
        array<int, 2> buf_sizes{M(bidx).target_shape()};
        buf_sizes() = buf_size;

        std::cout << "block_name = " << bname << " nfft_buf_size = " << buf_size << "\n";

        M_nfft(bidx) = nfft_array_t<2, 2>(M(bidx).mesh(), M(bidx).data(), buf_sizes);
      }
    }
  }

  template <G2_channel Channel> void measure_G2_iw_nfft<Channel>::accumulate(mc_weight_t s) {

    auto nfft_fill = [this](det_type const &det, nfft_array_t<2, 2> &nfft_matrix) {
      const double beta = this->data.config.beta();
      foreach (det, [&nfft_matrix, beta](op_t const &x, op_t const &y, det_scalar_t M) {
        nfft_matrix.push_back({beta - double(x.first), double(y.first)}, {x.second, y.second}, M);
      })
        ;
    };

    timer_M.start();
    // Intermediate M matrices for all blocks
    M() = 0;
    for (auto bidx : range(M.size())) {
      nfft_fill(data.dets[bidx], M_nfft(bidx));
      M_nfft(bidx).flush();
    }
    timer_M.stop();

    // ---------------------------------------------------------------
    // Recombine products of scattering matrices
    // to accumulate two particle quantities

    accumulate_G2(s);
  }

  template class measure_G2_iw_nfft<G2_channel::AllFermionic>;
  template class measure_G2_iw_nfft<G2_channel::PP>;
  template class measure_G2_iw_nfft<G2_channel::PH>;

} // namespace triqs_cthyb
