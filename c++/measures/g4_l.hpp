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
#include <triqs/experimental/nfft_array.hpp>

#include <triqs/utility/legendre.hpp>

#include "../qmc_data.hpp"

#include "util.hpp"

namespace cthyb {

  using namespace triqs::arrays;
  using namespace triqs::experimental;

  // Generates values of \tilde P_l(x(\tau_1-\tau_2))
  struct tilde_p_gen {
    triqs::utility::legendre_generator l_gen;
    double beta;
    double f;
    tilde_p_gen(double beta) : beta(beta) {}
    void reset(time_pt const &tau1, time_pt const &tau2) {
      l_gen.reset(2 * double(tau1 - tau2) / beta - 1);
      f = tau1 > tau2 ? 1 : -1;
    }
    double next() { return f * l_gen.next(); }
  };
  
  // Measure G^4(iomega,l,l')
  template <g4_channel Channel> struct measure_g4_l {

    qmc_data const &data;
    g4_wll_t::view_type g4_wll;
    g4_measures_t g4_measures;
    block_order order;

    mc_weight_t average_sign;

    // Object that performs NFFT transform
    array<nfft_array_t<1, 6>, 2> nfft_buf;

    measure_g4_l(std::optional<g4_wll_t> & g4_wll_opt, qmc_data const &data, g4_measures_t & g4_measures);
    void accumulate(mc_weight_t s);
    void collect_results(triqs::mpi::communicator const &c);

    // internal methods 
    double setup_times(tilde_p_gen & p_l1_gen, tilde_p_gen & p_l2_gen, op_t const & i, op_t const & j, op_t const & k, op_t const & l);
  };
}
