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
#include <triqs/gfs/types.hpp>
#include <triqs/utility/callbacks.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/statistics/histograms.hpp>

#include "container_set.hpp"
#include "solve_parameters.hpp"
#include "atom_diag.hpp"
#include "atom_diag_functions.hpp"

namespace cthyb {

  /// Core class of the cthyb solver
  class solver_core : public container_set_t {

    double beta;                                   // inverse temperature
    atom_diag h_diag;                              // diagonalization of the local problem
    gf_struct_t gf_struct; // Block structure of the Green function FIXME
    many_body_op_t _h_loc;                         // The local Hamiltonian = h_int + h0
    int n_iw, n_tau, n_l;

    // Single-particle Green's function containers
    g_iw_t _G0_iw; // Non-interacting Matsubara Green's function
    g_tau_t _Delta_tau; // Imaginary-time Hybridization function

    // Two-particle Green's functions
    //g4_tau_t _G2_tau; // Two-particle Green's function (three fermionic imaginary times)
    g4_iw_t _G2_inu; // Two-particle Green's function (three fermionic matsubaras)

    g4_iw_t _G2_iw_inu_inup_pp; // Two-particle Green's function, fermionic matsubaras, pp-channel
    g4_iw_t _G2_iw_inu_inup_ph; // Two-particle Green's function, fermionic matsubaras, ph-channel

    g4_wll_t _G2_iw_l_lp_pp; // Two-particle Green's function, Legendre coefficients, pp-channel
    g4_wll_t _G2_iw_l_lp_ph; // Two-particle Green's function, Legendre coefficients, ph-channel

    histogram _pert_order_total;               // Histogram of the total perturbation order
    histo_map_t _pert_order;                   // Histograms of the perturbation order for each block
    std::vector<matrix_t> _density_matrix;     // density matrix, when used in Norm mode
    triqs::mpi::communicator _comm;            // define the communicator, here MPI_COMM_WORLD
    solve_parameters_t _last_solve_parameters; // parameters of the last call to solve
    histo_map_t _performance_analysis;         // Histograms used for performance analysis
    mc_weight_t _average_sign;                 // average sign of the QMC
    int _solve_status;                         // Status of the solve upon exit: 0 for clean termination, > 0 otherwise.

    public:
    solver_core(double beta, std::map<std::string, indices_type> const &gf_struct, int n_iw = 1025, int n_tau = 10001, int n_l = 50);

    /// Solve the impurity problem for the given Hamiltonian h_loc and with specified parameters params.
    TRIQS_WRAP_ARG_AS_DICT // Wrap the solver parameters as a dictionary in python with the clang tool
       void
       solve(solve_parameters_t const &p);

    /// Return reference to the container_set
    container_set_t & container_set() { return static_cast<container_set_t &>(*this); }
    container_set_t const & container_set() const { return static_cast<container_set_t const &>(*this); }

    /// The local Hamiltonian of the problem: :math:`H_{loc}` used in the last call to ``solve()``.
    many_body_op_t const &h_loc() const { return _h_loc; }

    /// Set of parameters used in the last call to ``solve()``.
    solve_parameters_t last_solve_parameters() const { return _last_solve_parameters; }

    /// :math:`\Delta(\tau)` in imaginary time.
    block_gf_view<imtime> Delta_tau() { return _Delta_tau; }

    /// :math:`G_0(i\omega)` in imaginary frequencies.
    block_gf_view<imfreq> G0_iw() { return _G0_iw; }

    /// Accumulated :math:`G(\tau)` in imaginary time.
    //block_gf_view<imtime> G_tau() { return _G_tau; }

    /// Accumulated :math:`G_l` in Legendre polynomials representation.
    //block_gf_view<legendre> G_l() { return _G_l; }

    /// Accumulated two-particle Green's function :math:`G^{(2)}(\tau_1,\tau_2,\tau_3)`
    //block2_gf_view<cartesian_product<imtime, imtime, imtime>, tensor_valued<4>> G2_tau() { return _G2_tau; }

    /// Accumulated two-particle Green's function :math:`G^{(2)}(i\nu,i\nu',i\nu'')`
    //block2_gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> G2_inu() { return _G2_inu; }
    
    /// Accumulated two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the pp-channel.
    //block2_gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> G2_iw_inu_inup_pp() { return _G2_iw_inu_inup_pp; }

    /// Accumulated two-particle Green's function :math:`G^{(2)}(i\omega,i\nu,i\nu')` in the ph-channel.
    //block2_gf_view<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> G2_iw_inu_inup_ph() { return _G2_iw_inu_inup_ph; }

    /// Accumulated two-particle Green's function :math:`G^{(2)}(i\omega,l,l')` in the pp-channel.
    block2_gf_view<cartesian_product<imfreq, legendre, legendre>, tensor_valued<4>> G2_iw_l_lp_pp() { return _G2_iw_l_lp_pp; }

    /// Accumulated two-particle Green's function :math:`G^{(2)}(i\omega,l,l')` in the ph-channel.
    block2_gf_view<cartesian_product<imfreq, legendre, legendre>, tensor_valued<4>> G2_iw_l_lp_ph() { return _G2_iw_l_lp_ph; }

    /// Atomic :math:`G(\tau)` in imaginary time.
    block_gf_view<imtime> atomic_gf() const { return ::cthyb::atomic_gf(h_diag, beta, gf_struct, _Delta_tau[0].mesh().size()); }

    /// Accumulated density matrix.
    std::vector<matrix_t> const &density_matrix() const { return _density_matrix; }

    /// Diagonalization of :math:`H_{loc}`.
    atom_diag const &h_loc_diagonalization() const { return h_diag; }

    /// Histogram of the total perturbation order.
    histogram const &get_perturbation_order_total() const { return _pert_order_total; }

    /// Histograms of the perturbation order for each block.
    histo_map_t const &get_perturbation_order() const { return _pert_order; }

    /// Histograms related to the performance analysis.
    histo_map_t const &get_performance_analysis() const { return _performance_analysis; }

    /// Monte Carlo average sign.
    mc_weight_t average_sign() const { return _average_sign; }

    /// Status of the ``solve()`` on exit.
    int solve_status() const { return _solve_status; }
  };
}
