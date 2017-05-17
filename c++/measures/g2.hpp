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
#include <triqs/statistics/histograms.hpp>
#include "../nfft_array.hpp"
#include "../qmc_data.hpp"

#include "util.hpp"

namespace cthyb {

  using namespace triqs::arrays;

  // Measure one block of G^2(iomega,inu,inu')
  template <g2_channel Channel> struct measure_g2_inu {

    using g2_view_type = gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>>;
    using M_type = gf<cartesian_product<imfreq, imfreq>, matrix_valued>;
    using M_block_type = block_gf<cartesian_product<imfreq, imfreq>, matrix_valued>;

    qmc_data const &data;
    g2_view_type g2;
    const int A;           // Index of block A within gf_struct
    const int B;           // Index of block B within gf_struct
    const bool diag_block; // A == B
    mc_weight_t average_sign;
    block_order order;

    // Objects that perform 2D NFFT transforms
    nfft_array_t<2, 2> nfft_M_ab, nfft_M_cd;
    nfft_array_t<2, 2> nfft_M_ad, nfft_M_cb;

    // Results of NFFT transforms
    M_type M_ab, M_cd, M_ad, M_cb;

    M_block_type M;
    array<nfft_array_t<2, 2>, 1> M_nfft;

    measure_g2_inu(int A, int B, g2_view_type g2, qmc_data const &data, int buf_size_A, int buf_size_B, block_order order, gf_struct_t const & gf_struct);
    void accumulate(mc_weight_t s);
    void collect_results(triqs::mpi::communicator const &c);

    inline void accumulate_impl_AABB(mc_weight_t s, M_type const & M_ab, M_type const & M_cd);
    inline void accumulate_impl_ABBA(mc_weight_t s, M_type const & M_ad, M_type const & M_cb);

    // Frequency placeholders
    clef::placeholder<0> iw_;
    clef::placeholder<1> inu_;
    clef::placeholder<2> inup_;

    // Index placeholders
    clef::placeholder<3> a_;
    clef::placeholder<4> b_;
    clef::placeholder<5> c_;
    clef::placeholder<6> d_;

    clef::placeholder<0> w1;
    clef::placeholder<1> w2;
    clef::placeholder<2> w3;
  };

  // Measure one block of G^2(iomega,l,l')
  template <g2_channel Channel, block_order Order> struct measure_g2_legendre {

    using g2_view_type = gf_view<cartesian_product<imfreq, legendre, legendre>, tensor_valued<4>>;

    qmc_data const &data;
    g2_view_type g2;
    const int A;                 // Index of block A within gf_struct
    const int B;                 // Index of block B within gf_struct
    const bool diag_block;       // A == B
    const size_t n_l;            // Number of Legendre coefficients
    mc_weight_t z;
    int64_t num;

    // Object that performs NFFT transform
    nfft_array_t<1, 6> nfft_abcd;

    measure_g2_legendre(int b1, int b2, g2_view_type g2, qmc_data const &data, int buf_size_A, int buf_size_B);
    void accumulate(mc_weight_t s);
    void collect_results(triqs::mpi::communicator const &c);
  };
}
