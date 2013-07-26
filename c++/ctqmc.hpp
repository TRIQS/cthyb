/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by I. Krivenko, M. Ferrero, O. Parcollet
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
#ifndef TRIQS_CTQMC_KRYLOV_CTQMC_HPP
#define TRIQS_CTQMC_KRYLOV_CTQMC_HPP

#include <triqs/utility/first_include.hpp>
#include <triqs/parameters.hpp>
#include <triqs/gf/block.hpp>
#include <triqs/gf/imtime.hpp>
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>

#include "operator.hpp"
#include "qmc_data.hpp"

using triqs::utility::parameters;
using triqs::utility::parameter_defaults;
using triqs::utility::many_body_operator;
using triqs::mc_tools::mc_generic;
using triqs::gf::gf_view;
using triqs::gf::gf;
using triqs::gf::block_index;
using triqs::gf::imtime;

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

 class ctqmc {

  typedef std::complex<double> mc_sign_type;

  parameters apply_defaults(parameters const& p);
  
  public:

  typedef gf_view<block_index, gf::gf<imtime>> imtime_gf_t;

  template<typename ...IndexType>
  ctqmc (
    parameters const& p,
    many_body_operator<double, IndexType...> const& h_loc,
    std::vector<many_body_operator<double, IndexType...>> const & quantum_numbers,
    fundamental_operator_set<IndexType...> const & fops,
    std::map<std::tuple<IndexType...>, std::pair<int,int>> const & indices_to_ints,
    imtime_gf_t G_tau,
    imtime_gf_t Delta_tau
    ) :
  params(apply_defaults(p)),
  G_tau(G_tau),
  data(params,Delta_tau,sorted_spaces(h_loc,quantum_numbers,fops,indices_to_ints)),
  qmc(params) {}
  
  void solve();

  private:

  parameters params;
  imtime_gf_t G_tau;
  qmc_data data;
  mc_generic<mc_sign_type> qmc;
 };

}}}}

#endif
