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
#include <triqs/operators/many_body_operator.hpp>
#include "solve_parameters.hpp"
#include "./qmc_data.hpp"

namespace cthyb {

using namespace triqs::params;
using namespace triqs::utility;
using mc_weight_type = double;

class solver_core {

 using mc_sign_type = mc_weight_type;

 double beta;
 sorted_spaces sosp;
 std::map<std::string,std::vector<int>> gf_struct;
 block_gf<imfreq> _G0_iw;                  // Green's function containers: imaginary-freq Green's functions
 block_gf<imtime> _Delta_tau, _G_tau;      // Green's function containers: imaginary-time Green's functions
 block_gf<legendre> _G_l;                  // Green's function containers: Legendre coefficients
 boost::mpi::communicator _comm;           // define the communicator, here MPI_COMM_WORLD

 public:

 solver_core(double beta, std::map<std::string,std::vector<int>> const & gf_struct, int n_iw=1025, int n_tau=10001, int n_l=50);


 /// Solve the impurity problem for the given Hamiltonian h_loc and with specified parameters params.
 __attribute__((annotate("use_parameter_class")))
 void solve(solve_parameters_t const & p);

 // input containers
 /// G0(iw) in imaginary frequencies
 block_gf_view<imfreq> G0_iw() { return _G0_iw; }
 //void set_G0_iw(block_gf_view<imfreq> G0) { _G0_iw = G0; }

 /// Delta(tau) in imaginary time
 block_gf_view<imtime> Delta_tau() { return _Delta_tau; }

 // imaginary-time measurements
 /// G(tau) in imaginary time
 block_gf_view<imtime> G_tau() { return _G_tau; }
 
 // Legendre measurements
 /// G_l in Legendre polynomials representation
 block_gf_view<legendre> G_l() { return _G_l; }

 // Atomic GF (without hybridization)
 /// Atomic G(tau) in imaginary time
 block_gf_view<imtime> atomic_gf() const { return sosp.atomic_gf(beta,gf_struct,_G_tau[0].mesh().size()); }

};

// Change imaginary-time mesh of a GF
gf<imtime> change_mesh(gf_const_view<imtime> old_gf, int new_n_tau);
block_gf<imtime> change_mesh(block_gf_const_view<imtime> old_gf, int new_n_tau);

}
