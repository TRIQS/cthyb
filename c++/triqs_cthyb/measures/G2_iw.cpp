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

#include <boost/math/constants/constants.hpp>

#include "./G2_iw.hpp"

namespace triqs_cthyb {

  template <G2_channel Channel>
  measure_G2_iw<Channel>::measure_G2_iw(std::optional<G2_iw_t> &G2_iw_opt, qmc_data const &data, G2_measures_t const &G2_measures)
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

    // Allocate temporary two-frequency matrix M
    {
      int nfreq = std::max(3 * n_fermionic, n_bosonic + n_fermionic);
      gf_mesh<imfreq> iw_mesh_large{beta, Fermion, nfreq};
      //gf_mesh<cartesian_product<imfreq, imfreq>> M_mesh{iw_mesh_large, iw_mesh_large};
      M_mesh = M_mesh_t{iw_mesh_large, iw_mesh_large};

      if (Channel == G2_channel::AllFermionic) { // Smaller mesh possible in AllFermionic
        gf_mesh<imfreq> iw_mesh_small{beta, Fermion, n_fermionic};
        M_mesh = gf_mesh<cartesian_product<imfreq, imfreq>>{iw_mesh_large, iw_mesh_small};
      }

      // Initialize intermediate scattering matrix
      M = block_gf{M_mesh, G2_measures.gf_struct};

      for (auto const &m : M) {
        auto norb1        = static_cast<size_t>(m.target_shape()[0]);
        auto norb2        = static_cast<size_t>(m.target_shape()[1]);
        size_t nfreq_pts = static_cast<size_t>(std::get<0>(M_mesh.components()).full_size());
        M_block_arr.push_back(array<std::complex<double>, 4>{norb1, norb2, nfreq_pts, nfreq_pts});
      }
    }
  }

  namespace {

    // Index placeholders
    clef::placeholder<0> i;
    clef::placeholder<1> j;
    clef::placeholder<2> k;
    clef::placeholder<3> l;

    // Frequency placeholders
    clef::placeholder<4> w;
    clef::placeholder<5> n1;
    clef::placeholder<6> n2;
    clef::placeholder<7> n3;

  } // namespace

  template <G2_channel Channel> void measure_G2_iw<Channel>::accumulate(mc_weight_t s) {

    s *= data.atomic_reweighting;
    average_sign += s;

    auto M_ww_fill = [this](det_type const &det, M_t &M_ww) {
      const double beta = this->data.config.beta();
      foreach (det, [&M_ww, beta](op_t const &x, op_t const &y, det_scalar_t M_xy) {
        // insert accumulation
        double t1 = double(x.first);
        double t2 = double(y.first);

        for (auto const & [ w1, w2 ] : M_ww.mesh()) {
	  M_ww[w1, w2](x.second, y.second) += exp((beta - t1) * w1) * M_xy * exp(t2 * w2);
	}
      })
        ;
    };

    {
    timer_M_ww_fill.start();
    // Intermediate M matrices for all blocks
    M() = 0;
    for (auto bidx : range(M.size())) { M_ww_fill(data.dets[bidx], M[bidx]); }
    timer_M_ww_fill.stop();
    }
    
    // ========================================================

    const double beta = data.config.beta();
    const double pi_beta = boost::math::constants::pi<double>() / beta;
    
    auto M_arr_fill = [pi_beta, beta](det_type const &det, M_arr_t &M_arr, M_mesh_t const & M_mesh) {
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
	
	for (auto const i1 : range(M_arr.shape()[2]) ) {

	  auto exp2 = exp(dWt2 * (mesh2.first_index() + 0.5));
	  auto exp1_M_xy = exp1 * M_xy;

	  for (auto const i2 : range(M_arr.shape()[3]) ) {

	    M_arr(x.second, y.second, i1, i2) += exp1_M_xy * exp2;
	    exp2 *= dexp2;
	    
	  }
	  exp1 *= dexp1;
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

    for (auto bidx : range(M_block_arr.size())) {
      M[bidx].data()(n1, n2, i, j) << M_block_arr[bidx](i, j, n1, n2);

      /*
      auto M_arr = M_block_arr[bidx];

      for (auto const a : range(M_arr.shape()[0]) )
      for (auto const b : range(M_arr.shape()[1]) )
      for (auto const w1 : range(M_arr.shape()[2]) )
      for (auto const w2 : range(M_arr.shape()[3]) )
	M[bidx].data()(w1, w2, a, b) = M_block_arr[bidx](a, b, w1, w2);
      //M[bidx][{w1, w2}](a, b) = M_block_arr[bidx](a, b, w1, w2);
      */
    }

    timer_M_arr_fill.stop();
    }

    // ========================================================
    {
      timer_MM_prod.start();
    for (auto &m : G2_measures()) {
      auto G2_iw_block = G2_iw(m.b1.idx, m.b2.idx);
      bool diag_block  = (m.b1.idx == m.b2.idx);
      if (order == block_order::AABB || diag_block) accumulate_impl_AABB(G2_iw_block, s, M(m.b1.idx), M(m.b2.idx));
      if (order == block_order::ABBA || diag_block) accumulate_impl_ABBA(G2_iw_block, s, M(m.b1.idx), M(m.b2.idx));
    }
      timer_MM_prod.stop();
    }
    
  }

  // -- Particle-hole

  template <>
  inline void measure_G2_iw<G2_channel::PH>::accumulate_impl_AABB(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) {
    G2(w, n1, n2)
    (i, j, k, l) << G2(w, n1, n2)(i, j, k, l)                    //
          + s * M_ij(n1, n1 + w)(i, j) * M_kl(n2 + w, n2)(k, l); // sign in lhs in fft
  }

  template <>
  inline void measure_G2_iw<G2_channel::PH>::accumulate_impl_ABBA(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_il, M_t const &M_kj) {
    G2(w, n1, n2)
    (i, j, k, l) << G2(w, n1, n2)(i, j, k, l)                    //
          - s * M_il(n1, n2)(i, l) * M_kj(n2 + w, n1 + w)(k, j); // sign in lhs in fft
  }

  // -- Particle-particle

  template <>
  inline void measure_G2_iw<G2_channel::PP>::accumulate_impl_AABB(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) {
    G2(w, n1, n2)
    (i, j, k, l) << G2(w, n1, n2)(i, j, k, l)                    //
          + s * M_ij(n1, w - n2)(i, j) * M_kl(w - n1, n2)(k, l); // sign in lhs in fft
  }

  template <>
  inline void measure_G2_iw<G2_channel::PP>::accumulate_impl_ABBA(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_il, M_t const &M_kj) {
    G2(w, n1, n2)
    (i, j, k, l) << G2(w, n1, n2)(i, j, k, l)                    //
          - s * M_il(n1, n2)(i, l) * M_kj(w - n1, w - n2)(k, j); // sign in lhs in fft
  }

  // -- Fermionic

  template <>
  inline void measure_G2_iw<G2_channel::AllFermionic>::accumulate_impl_AABB(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij,
                                                                            M_t const &M_kl) {

    int size_ij = M_ij.target_shape()[0];
    int size_kl = M_kl.target_shape()[0];

    auto const iw_mesh = std::get<0>(G2.mesh());
    using mesh_point_t = typename decltype(iw_mesh)::mesh_point_t;

    for (auto const &n1 : iw_mesh)
      for (auto const &n2 : iw_mesh)
        for (auto const &n3 : iw_mesh) {
          mesh_point_t n4{iw_mesh, n1.index() + n3.index() - n2.index()};
          for (int i : range(size_ij))
            for (int j : range(size_ij))
              for (int k : range(size_kl))
                for (int l : range(size_kl)) G2[{n1, n2, n3}](i, j, k, l) += s * M_ij[{n2, n1}](j, i) * M_kl[{n4, n3}](l, k);
        }

    //G2(n1, n2, n3)(i, j, k, l) << G2(n1, n2, n3)(i, j, k, l) + s * M_ij(n2, n1)(j, i) * M_kl(n1 + n3 - n2, n3)(l, k);
  }

  template <>
  inline void measure_G2_iw<G2_channel::AllFermionic>::accumulate_impl_ABBA(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_il,
                                                                            M_t const &M_kj) {

    int size_il = M_il.target_shape()[0];
    int size_kj = M_kj.target_shape()[0];

    auto const iw_mesh = std::get<0>(G2.mesh());
    using mesh_point_t = typename decltype(iw_mesh)::mesh_point_t;

    for (auto const &n1 : iw_mesh)
      for (auto const &n2 : iw_mesh)
        for (auto const &n3 : iw_mesh) {
          mesh_point_t n4{iw_mesh, n1.index() + n3.index() - n2.index()};
          for (int i : range(size_il))
            for (int j : range(size_kj))
              for (int k : range(size_kj))
                for (int l : range(size_il)) G2[{n1, n2, n3}](i, j, k, l) += -s * M_il[{n4, n1}](l, i) * M_kj[{n2, n3}](j, k);
        }

    //G2(n1, n2, n3)(i, j, k, l) << G2(n1, n2, n3)(i, j, k, l) - s * M_il(n1 + n3 - n2, n1)(l, i) * M_kj(n2, n3)(j, k);
  }

  // --

  template <G2_channel Channel> void measure_G2_iw<Channel>::collect_results(triqs::mpi::communicator const &com) {
    average_sign = mpi_all_reduce(average_sign, com);
    G2_iw        = mpi_all_reduce(G2_iw, com);
    G2_iw = G2_iw / (real(average_sign) * data.config.beta());

    std::cout << "timer_M_arr_fill = " << double(timer_M_arr_fill) << "\n";    
    std::cout << "timer_M_ww_fill = " << double(timer_M_ww_fill) << "\n";    
    std::cout << "timer_MM_prod = " << double(timer_MM_prod) << "\n";    
  }

  template class measure_G2_iw<G2_channel::AllFermionic>;
  template class measure_G2_iw<G2_channel::PP>;
  template class measure_G2_iw<G2_channel::PH>;

} // namespace triqs_cthyb
