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
#include "qmc_data.hpp"

namespace cthyb {

using namespace triqs::arrays;

enum g2_channel {PP, PH};

namespace details {

// NFFT transform of a matrix-valued function of two tau arguments
class nfft_matrix_t {

// Perform NFFT checks
#ifdef NDEBUG
 static constexpr bool do_checks = false;
#else
 static constexpr bool do_checks = true;
#endif

public:

 using res_gf_t = gf<cartesian_product<imfreq,imfreq>, scalar_valued>;
 using res_gf_view = gf_const_view<cartesian_product<imfreq,imfreq>, scalar_valued>;

 nfft_matrix_t() = default;

 // size1, size2 - matrix dimensions
 // beta - inverse temperature
 // n_freq1, n_freq2 - number of positive Matsubara frequencies
 // to calculate the Fourier transformation on
 nfft_matrix_t(int size1, int size2, double beta, int n_freq1, int n_freq2);

 // Clear all input buffers
 void reset();

 // Add a new matrix element to the input buffer
 void push_back(std::pair<time_pt,int> const& x, std::pair<time_pt,int> const& y, dcomplex fxy);

 // Run transformation
 void transform();

 TRIQS_CLEF_IMPLEMENT_LAZY_CALL();
 // Access a matrix element of the transformation result
 res_gf_view operator()(int n1, int n2) const;

private:

 // Inverse temperature
 double beta;
 // Input and output data
 struct data_t {
  std::vector<std::pair<std::array<double,2>,dcomplex>> input;
  gf<cartesian_product<imfreq,imfreq>, scalar_valued> result;
 };
 array<data_t,2> data;
 // Maximal number of tau-pairs over all matrix elements
 size_t max_n_tau;
};

}

// Measure one block of G^2(iomega,inu,inu')
template<g2_channel Channel, block_order Order> struct measure_g2_inu {

 qmc_data const& data;
 g2_iw_inu_inup_block_view g2;
 const int A;                 // Index of block A within gf_struct
 const int B;                 // Index of block B within gf_struct
 const bool diag_block;       // A == B
 mc_weight_t z;
 int64_t num;

 // Objects that perform 2D NFFT transforms
 details::nfft_matrix_t nfft_matrix_ab, nfft_matrix_cd;
 details::nfft_matrix_t nfft_matrix_ad, nfft_matrix_cb;

 measure_g2_inu(int A, int B, g2_iw_inu_inup_block_view g2, qmc_data const& data);
 void accumulate(mc_weight_t s);
 void collect_results(triqs::mpi::communicator const& c);

};

// Measure one block of G^2(iomega,l,l')
template<g2_channel Channel, block_order Order> struct measure_g2_legendre {

 qmc_data const& data;
 g2_iw_l_lp_block_view g2;
 const int A;                 // Index of block A within gf_struct
 const int B;                 // Index of block B within gf_struct
 const bool diag_block;       // A == B
 const size_t size_A, size_B; // Sizes of blocks A and B
 mc_weight_t z;
 int64_t num;

 measure_g2_legendre(int b1, int b2, g2_iw_l_lp_block_view g2, qmc_data const& data);
 void accumulate(mc_weight_t s);
 void collect_results(triqs::mpi::communicator const& c);

};

}
