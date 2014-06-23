/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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
#include <triqs/mc_tools.hpp>
#include <triqs/utility/callbacks.hpp>
#include <boost/mpi/communicator.hpp>

#include "./qmc_data.hpp"

namespace cthyb {

using namespace triqs::params;
using namespace triqs::utility;
using mc_weight_type = double;
using parameters_t = triqs::params::parameters;

class ctqmc {

 using mc_sign_type = mc_weight_type;

 double beta;
 sorted_spaces sosp;
 std::map<std::string,std::vector<int>> gf_struct;
 block_gf<imtime> deltat, gt;             // Green's function containers: imaginary-time Green's functions
 block_gf<imfreq> g0w;                    // Green's function containers: imaginary-freq Green's functions
 boost::mpi::communicator _comm;          // define the communicator, here MPI_COMM_WORLD

 public:
 using real_operator_t = many_body_operator<double>;

 ctqmc(double beta_, std::map<std::string,std::vector<int>> const & gf_struct, int n_iw=1025, int n_tau=10001);

 void solve(real_operator_t h_loc, params::parameters params,
            std::vector<real_operator_t> const & quantum_numbers = std::vector<real_operator_t> {},
            bool use_quantum_numbers = false);

 // input containers
 block_gf_view<imfreq> g0w_view() { return g0w; }
 block_gf_view<imtime> deltat_view() { return deltat; }

 // imaginary-time measurements
 block_gf_view<imtime> gt_view() { return gt; }

 // specify all required and optional parameters and generate help from them
 static parameters solve_parameters();
 static void help();

 block_gf_view<imtime> atomic_gf(double beta_) const { return sosp.atomic_gf(beta_); }

};
}
