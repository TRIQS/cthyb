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
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/statistics/histograms.hpp>

#include "types.hpp"
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

    /// The local Hamiltonian of the problem: :math:`H_{loc}` used in the last call to ``solve()``.
    many_body_op_t const &h_loc() const { return _h_loc; }

    /// Set of parameters used in the last call to ``solve()``.
    solve_parameters_t last_solve_parameters() const { return _last_solve_parameters; }

    /// :math:`\Delta(\tau)` in imaginary time.
    block_gf_view<imtime> Delta_tau() { return _Delta_tau; }

    /// :math:`G_0(i\omega)` in imaginary frequencies.
    block_gf_view<imfreq> G0_iw() { return _G0_iw; }

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
