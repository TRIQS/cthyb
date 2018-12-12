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
#include "./G2_iw_acc.hpp"

namespace triqs_cthyb {

  template <G2_channel Channel>
  measure_G2_iw_nfft<Channel>::measure_G2_iw_nfft(std::optional<G2_iw_t> &G2_iw_opt, qmc_data const &data, G2_measures_t const &G2_measures)
     : data(data), average_sign(0), G2_measures(G2_measures) {

    const double beta = data.config.beta();

    order           = G2_measures.params.measure_G2_block_order;
    int n_bosonic   = G2_measures.params.measure_G2_n_bosonic;
    int n_fermionic = G2_measures.params.measure_G2_n_fermionic;

    // Allocate the two-particle Green's function
    {
      gf_mesh<imfreq> mesh_f{beta, Fermion, n_fermionic};
      gf_mesh<imfreq> mesh_b{beta, Boson, n_bosonic};

      gf_mesh<cartesian_product<imfreq, imfreq, imfreq>> mesh_fff{mesh_f, mesh_f, mesh_f};
      gf_mesh<cartesian_product<imfreq, imfreq, imfreq>> mesh_bff{mesh_b, mesh_f, mesh_f};

      if (Channel == G2_channel::AllFermionic)
        G2_iw_opt = make_block2_gf(mesh_fff, G2_measures.gf_struct, order);
      else
        G2_iw_opt = make_block2_gf(mesh_bff, G2_measures.gf_struct, order);

      G2_iw.rebind(*G2_iw_opt);
      G2_iw() = 0;
    }

    // Allocate temporary NFFT two-frequency matrix M
    {
      gf_mesh<cartesian_product<imfreq, imfreq>> M_mesh;

      if (Channel == G2_channel::AllFermionic) { // Smaller mesh possible in AllFermionic
	gf_mesh<imfreq> iw_mesh_large{beta, Fermion, 3 * n_fermionic};
        gf_mesh<imfreq> iw_mesh_small{beta, Fermion, n_fermionic};
        M_mesh = gf_mesh<cartesian_product<imfreq, imfreq>>{iw_mesh_large, iw_mesh_small};
      } else {
	int nfreq = n_bosonic + n_fermionic;
	gf_mesh<imfreq> iw_mesh{beta, Fermion, nfreq};
        M_mesh = gf_mesh<cartesian_product<imfreq, imfreq>>{iw_mesh, iw_mesh};
      }      

      // Initialize intermediate scattering matrix
      M = block_gf(M_mesh, G2_measures.gf_struct);
    }

    // Initialize the nfft_buffers mirroring the matrix M
    {
      M_nfft.resize(M.size());

      for (auto bidx : range(M.size())) {
        std::string bname = M.block_names()[bidx];

        int buf_size = 10; // default
        if (G2_measures.params.nfft_buf_sizes.count(bname)) { buf_size = G2_measures.params.nfft_buf_sizes.at(bname); }
        array<int, 2> buf_sizes{M(bidx).target_shape()};
        buf_sizes() = buf_size;

	std::cout << "bname = " << bname << " buf_size = " << buf_size << "\n";
	
        M_nfft(bidx) = nfft_array_t<2, 2>(M(bidx).mesh(), M(bidx).data(), buf_sizes);
      }
    }
  }

  template <G2_channel Channel> void measure_G2_iw_nfft<Channel>::accumulate(mc_weight_t s) {

    s *= data.atomic_reweighting;
    average_sign += s;

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
    timer_G2.start();

    for (auto &m : G2_measures()) {
      auto G2_iw_block = G2_iw(m.b1.idx, m.b2.idx);
      bool diag_block  = (m.b1.idx == m.b2.idx);

      /*
      if(Channel == G2_channel::PH && m.target_shape[0] == 1) {
	if (order == block_order::AABB || diag_block) accumulate_impl_AABB_opt<Channel>(G2_iw_block, s, M(m.b1.idx), M(m.b2.idx));
      if (order == block_order::ABBA || diag_block) accumulate_impl_ABBA_opt<Channel>(G2_iw_block, s, M(m.b1.idx), M(m.b2.idx));
      } else {
      */

      if (order == block_order::AABB || diag_block) accumulate_impl_AABB<Channel>(G2_iw_block, s, M(m.b1.idx), M(m.b2.idx));
      if (order == block_order::ABBA || diag_block) accumulate_impl_ABBA<Channel>(G2_iw_block, s, M(m.b1.idx), M(m.b2.idx));

      //}
      
    }
    
    timer_G2.stop();
  }

  template <G2_channel Channel> void measure_G2_iw_nfft<Channel>::collect_results(triqs::mpi::communicator const &com) {
    average_sign = mpi_all_reduce(average_sign, com);
    G2_iw        = mpi_all_reduce(G2_iw, com);
    
    for (auto &g2_iw : G2_iw) g2_iw /= (real(average_sign) * data.config.beta());
    //G2_iw = G2_iw / (real(average_sign) * data.config.beta()); // This segfaults on triqs/unstable da793fbd

    std::cout << "measure/G2_iw_nfft: timer_M  = " << double(timer_M) << "\n";
    std::cout << "measure/G2_iw_nfft: timer_G2 = " << double(timer_G2) << "\n";
  }

  template class measure_G2_iw_nfft<G2_channel::AllFermionic>;
  template class measure_G2_iw_nfft<G2_channel::PP>;
  template class measure_G2_iw_nfft<G2_channel::PH>;

} // namespace triqs_cthyb
