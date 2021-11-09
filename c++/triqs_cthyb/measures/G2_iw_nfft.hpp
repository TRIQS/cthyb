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
#pragma once

#include <triqs/experimental/nfft_array.hpp>

#include "G2_iw_acc.hpp"

namespace triqs_cthyb {

  using namespace nda;
  using namespace triqs::experimental;

  // Measure the two-particle Green's function in Matsubara frequency
  template <G2_channel Channel> class measure_G2_iw_nfft : public G2_iw::measure_G2_iw_base<Channel> {

    public:
    measure_G2_iw_nfft(std::optional<G2_iw_t> &G2_iw_opt, qmc_data const &data, G2_measures_t const &G2_measures);
    void accumulate(mc_weight_t s);

    using B = G2_iw::measure_G2_iw_base<Channel>;
    using B::collect_results;
    
    private:
    array<nfft_array_t<2, 2>, 1> M_nfft;
    using B::M, B::M_mesh, B::G2_measures, B::data, B::timer_M, B::accumulate_G2;
  };

} // namespace triqs_cthyb
