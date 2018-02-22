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

#include "./G2_iw_w0.hpp"
#include "./G2_iw_acc.hpp"

#include <boost/math/constants/constants.hpp>

namespace triqs_cthyb {

  template <G2_channel Channel>
  measure_G2_iw_w0<Channel>::measure_G2_iw_w0(std::optional<G2_iw_t> &G2_iw_opt, qmc_data const &data, G2_measures_t const &G2_measures)
     : data(data), average_sign(0), G2_measures(G2_measures) {
    
    TRIQS_RUNTIME_ERROR << "G2_iw_w0 measure not working yet.";

     //static_assert(Channel != G2_channel::AllFermionic,
    //	  "G2_iw_w0 measure, zero bosonic transfer not possible for G2_Channel::AllFermionic");

    if(Channel == G2_channel::AllFermionic)
      TRIQS_RUNTIME_ERROR << "G2_iw_w0 measure, zero bosonic transfer not possible for G2_Channel::AllFermionic";
      
    const double beta = data.config.beta();

    order           = G2_measures.params.measure_G2_block_order;
    int n_fermionic = G2_measures.params.measure_G2_n_fermionic;

    // Allocate the two-particle Green's function
    {
      gf_mesh<imfreq> mesh_f{beta, Fermion, n_fermionic};
      gf_mesh<imfreq> mesh_b{beta, Boson, 1};

      gf_mesh<cartesian_product<imfreq, imfreq, imfreq>> mesh_bff{mesh_b, mesh_f, mesh_f};

      G2_iw_opt = make_block2_gf(mesh_bff, G2_measures.gf_struct, order);

      G2_iw.rebind(*G2_iw_opt);
      G2_iw() = 0;
    }

    // Allocate temporary two-frequency matrix M

    /*
    {
      int nfreq = n_bosonic + n_fermionic;
      M_mesh    = gf_mesh<imfreq>{beta, Fermion, n_fermionic};
    }
    */

      if (Channel == G2_channel::AllFermionic) { // Smaller mesh possible in AllFermionic
	gf_mesh<imfreq> iw_mesh_large{beta, Fermion, 3 * n_fermionic};
        gf_mesh<imfreq> iw_mesh_small{beta, Fermion, n_fermionic};
        M_mesh = gf_mesh<cartesian_product<imfreq, imfreq>>{iw_mesh_large, iw_mesh_small};
      } else {
	int n_bosonic = 1;
	int nfreq = n_bosonic + n_fermionic;
	gf_mesh<imfreq> iw_mesh{beta, Fermion, nfreq};
        M_mesh = gf_mesh<cartesian_product<imfreq, imfreq>>{iw_mesh, iw_mesh};
      }    

    // Initialize intermediate scattering matrix
    M = make_block_gf(M_mesh, G2_measures.gf_struct);

    // Accumulation buffer for scattering matrix
    for (auto const &m : M) {
      auto norb1       = static_cast<size_t>(m.target_shape()[0]);
      auto norb2       = static_cast<size_t>(m.target_shape()[1]);
      size_t nfreq_pts = static_cast<size_t>(std::get<0>(M_mesh).full_size());
      M_block_arr.push_back(M_diag_arr_t{norb1, norb2, nfreq_pts});
    }
  }

  template <G2_channel Channel> void measure_G2_iw_w0<Channel>::accumulate(mc_weight_t s) {

    s *= data.atomic_reweighting;
    average_sign += s;

    const double beta    = data.config.beta();
    const double pi_beta = boost::math::constants::pi<double>() / beta;

    // -- specialization used for zero bosonic transfer in PH-channel

    auto M_arr_fill = [pi_beta, beta](det_type const &det, M_diag_arr_t &M_arr, M_mesh_t const &M_mesh) {
      foreach (det, [&M_mesh, &M_arr, pi_beta, beta](op_t const &x, op_t const &y, det_scalar_t M_xy) {

        double t1 = double(x.first);
        double t2 = double(y.first);

        auto mesh1 = std::get<0>(M_mesh.components());
        auto mesh2 = std::get<1>(M_mesh.components());

        std::complex<double> dWt1(0., 2 * pi_beta * (beta - t1));
        std::complex<double> dWt2(0., 2 * pi_beta * t2);

        auto dexp1 = exp(dWt1);
        auto dexp2 = exp(dWt2);

        auto exp1 = exp(dWt1 * (mesh1.first_index() + 0.5));
        auto exp2 = exp(dWt2 * (mesh2.first_index() + 0.5));

        for (auto const i1 : range(M_arr.shape()[2])) {
          M_arr(x.second, y.second, i1) += exp1 * M_xy * exp2;
          exp1 *= dexp1;
          exp2 *= dexp2;
        }

      })
        ;
    };

    {
      timer_M_arr_fill.start();
      // Intermediate M matrices for all blocks
      for (auto bidx : range(M_block_arr.size())) {
        M_block_arr[bidx] *= 0;
        M_arr_fill(data.dets[bidx], M_block_arr[bidx], M_mesh);
      }

      // Reshuffle the accumulated scattering matrix into a Green's function object
      for (auto bidx : range(M_block_arr.size())) {
        {
          clef::placeholder<0> n1;
          clef::placeholder<1> n2;
          clef::placeholder<2> i;
          clef::placeholder<3> j;

          M[bidx].data()(n1, n2, i, j) << kronecker(n1, n2) * M_block_arr[bidx](i, j, n1);
        }
      }

      timer_M_arr_fill.stop();
    }

    // ---------------------------------------------------------------
    // Recombine products of scattering matrices
    // to accumulate two particle quantities

    {
      timer_MM_prod.start();

      for (auto &m : G2_measures()) {
        auto G2_iw_block = G2_iw(m.b1.idx, m.b2.idx);
        bool diag_block  = (m.b1.idx == m.b2.idx);
        if (Channel == G2_channel::PH && m.target_shape[0] == 1) {
          if (order == block_order::AABB || diag_block) accumulate_impl_AABB_opt<Channel>(G2_iw_block, s, M(m.b1.idx), M(m.b2.idx));
          if (order == block_order::ABBA || diag_block) accumulate_impl_ABBA_opt<Channel>(G2_iw_block, s, M(m.b1.idx), M(m.b2.idx));
        } else {
          if (order == block_order::AABB || diag_block) accumulate_impl_AABB<Channel>(G2_iw_block, s, M(m.b1.idx), M(m.b2.idx));
          if (order == block_order::ABBA || diag_block) accumulate_impl_ABBA<Channel>(G2_iw_block, s, M(m.b1.idx), M(m.b2.idx));
        }
      }

      timer_MM_prod.stop();
    }
  }

  template <G2_channel Channel> void measure_G2_iw_w0<Channel>::collect_results(triqs::mpi::communicator const &com) {
    average_sign = mpi_all_reduce(average_sign, com);
    G2_iw        = mpi_all_reduce(G2_iw, com);
    for (auto &g2_iw : G2_iw) g2_iw /= (real(average_sign) * data.config.beta());
    // G2_iw = G2_iw / (real(average_sign) * data.config.beta()); // This segfaults on triqs/unstable da793fbd

    std::cout << "G2_iw_w0 timer_M_arr_fill = " << double(timer_M_arr_fill) << "\n";
    std::cout << "G2_iw_w0 timer_M_ww_fill = " << double(timer_M_ww_fill) << "\n";
    std::cout << "G2_iw_w0 timer_MM_prod = " << double(timer_MM_prod) << "\n";
  }

  template class measure_G2_iw_w0<G2_channel::AllFermionic>;
  template class measure_G2_iw_w0<G2_channel::PP>;
  template class measure_G2_iw_w0<G2_channel::PH>;

} // namespace triqs_cthyb
