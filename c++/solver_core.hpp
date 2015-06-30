/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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

using namespace triqs::utility;
using mc_weight_type = double;
using mc_sign_type = mc_weight_type;
using indices_type = many_body_operator<double>::indices_t;
using many_body_op_type = triqs::utility::many_body_operator<double>;

class solver_core {

 double beta;
 sorted_spaces sosp;
 std::map<std::string, indices_type> gf_struct;
 many_body_op_type _h_loc;                  // The local Hamiltonian = h_int + h0
 block_gf<imfreq> _G0_iw;                   // Green's function containers: imaginary-freq Green's functions
 block_gf<imtime> _Delta_tau, _G_tau;       // Green's function containers: imaginary-time Green's functions
 block_gf<legendre> _G_l;                   // Green's function containers: Legendre coefficients
 boost::mpi::communicator _comm;            // define the communicator, here MPI_COMM_WORLD
 solve_parameters_t _last_solve_parameters; // parameters of the last call to solve
 mc_sign_type _average_sign;
 arrays::vector<double> state_trace_contribs;

 public:

 solver_core(double beta, std::map<std::string, indices_type> const & gf_struct, int n_iw=1025, int n_tau=10001, int n_l=50);

 /// Solve the impurity problem for the given Hamiltonian h_loc and with specified parameters params.
 TRIQS_WRAP_ARG_AS_DICT // Wrap the solver parameters as a dictionary in python with the clang tool
 void solve(solve_parameters_t const & p);

 /// The local Hamiltonian of the problem
 many_body_op_type h_loc() const { return _h_loc; }

 /// Set of parameters used in the last call to solve
 solve_parameters_t get_last_solve_parameters() const {return _last_solve_parameters;}

 /// G0(iw) in imaginary frequencies
 block_gf_view<imfreq> G0_iw() { return _G0_iw; }

 /// Delta(tau) in imaginary time
 block_gf_view<imtime> Delta_tau() { return _Delta_tau; }

 /// G(tau) in imaginary time
 block_gf_view<imtime> G_tau() { return _G_tau; }
 
 /// G_l in Legendre polynomials representation
 block_gf_view<legendre> G_l() { return _G_l; }

 /// Atomic G(tau) in imaginary time
 block_gf_view<imtime> atomic_gf() const { return sosp.atomic_gf(beta,gf_struct,_G_tau[0].mesh().size()); }

 /// Average contribution of each atomic state to the trace
 arrays::vector<double> const& get_state_trace_contribs() const { return state_trace_contribs; }

 /// Static observables of the atomic problem
 std::map<std::string,std::vector<double>> atomic_observables(std::map<std::string,real_operator_t> const& obs_map) const {
  return sosp.atomic_observables(obs_map);
 }

 /// Monte Carlo average sign
 mc_sign_type average_sign() const { return _average_sign; }

 /// Eigensystems of the atomic problem
 /// Returns a list of pairs (E,U), where H = U * diag(E) * U^+ (for each subspace)
 std::vector<std::pair<vector<double>,matrix<double>>> get_eigensystems() const {
  decltype(get_eigensystems()) res;
  for(auto const& es : sosp.get_eigensystems()) res.emplace_back(es.eigenvalues,es.unitary_matrix);
  return res;
 }

};

}
