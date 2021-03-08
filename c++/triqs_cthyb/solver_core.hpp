
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014-2017, H. U.R. Strand, P. Seth, I. Krivenko, 
 *                          M. Ferrero and O. Parcollet
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
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/atom_diag/functions.hpp>
#include <optional>

#include "types.hpp"
#include "container_set.hpp"
#include "parameters.hpp"

namespace triqs_cthyb {

  /// Core class of the cthyb solver
  class solver_core : public container_set_t {

    double beta;           // inverse temperature
    atom_diag h_diag;      // diagonalization of the local problem
    gf_struct_t gf_struct; // Block structure of the Green function
    many_body_op_t _h_loc; // The local Hamiltonian = h_int + h0
    int n_iw, n_tau, n_l;
    bool Delta_interface;

    std::vector<matrix_t> _density_matrix; // density matrix, when used in Norm mode
    mpi::communicator _comm;               // define the communicator, here MPI_COMM_WORLD
    histo_map_t _performance_analysis;     // Histograms used for performance analysis
    mc_weight_t _average_sign;             // average sign of the QMC
    int _solve_status;                     // Status of the solve upon exit: 0 for clean termination, > 0 otherwise.

    // Single-particle Green's function containers
    std::optional<G_iw_t> _G0_iw;      // Non-interacting Matsubara Green's function
    G_tau_t _Delta_tau; // Imaginary-time Hybridization function
    std::vector<matrix<dcomplex>> Delta_infty_vec; // Quadratic instantaneous part of G0_iw

    // Return reference to container_set
    container_set_t &container_set() { return static_cast<container_set_t &>(*this); }
    container_set_t const &container_set() const { return static_cast<container_set_t const &>(*this); }
 
    public:

    // Struct containing the parameters relevant for the solver construction
    constr_parameters_t constr_parameters;

    // Struct containing the parameters of the last call to the solve method
    solve_parameters_t solve_parameters;

    /**
     * Construct a CTHYB solver
     *
     * @param p Set of parameters specific to the CTHYB solver
     */
    CPP2PY_ARG_AS_DICT
    solver_core(constr_parameters_t const &p);

    // Delete assignement operator because of const members
    solver_core(solver_core const &p) = default;
    solver_core(solver_core &&p)      = default;
    solver_core &operator=(solver_core const &p) = delete;
    solver_core &operator=(solver_core &&p) = default;

    /**
     * Solve method that performs CTHYB calculation
     *
     * @param p Set of parameters for the CTHYB calculation
     */
    CPP2PY_ARG_AS_DICT
    void solve(solve_parameters_t const &p);

    /// The local Hamiltonian of the problem: :math:`H_{loc}` used in the last call to ``solve()``.
    many_body_op_t const &h_loc() const { return _h_loc; }

    /// Set of parameters used in the construction of the ``solver_core`` class.
    constr_parameters_t last_constr_parameters() const { return constr_parameters; }

    /// Set of parameters used in the last call to ``solve()``.
    solve_parameters_t last_solve_parameters() const { return solve_parameters; }

    /// :math:`G_0^{-1}(i\omega_n = \infty)` in Matsubara Frequency.
    std::vector<matrix<dcomplex>> Delta_infty() { return Delta_infty_vec; }

    /// Get a copy of the last container set.
    // HACK TO GET CPP2PY TO WRAP THE container_set_t struct.
    /*
    CPP2PY_ARG_AS_DICT
    void set_container_set(container_set_t &cs) { static_cast<container_set_t &>(*this) = cs; }
    container_set_t last_container_set() { return static_cast<container_set_t>(*this); }
    */
    
    /// :math:`\Delta(\tau)` in imaginary time.
    block_gf_view<imtime> Delta_tau() { return _Delta_tau; }

    /// :math:`G_0(i\omega)` in imaginary frequencies.
    block_gf_view<imfreq> G0_iw() { return _G0_iw.value(); }

    /// Atomic :math:`G(\tau)` in imaginary time.
    //block_gf_view<imtime> atomic_gf() const { return ::triqs_cthyb::atomic_gf(h_diag, beta, gf_struct, _Delta_tau[0].mesh().size()); }

    /// Accumulated density matrix.
    std::vector<matrix_t> const &density_matrix() const { return _density_matrix; }

    /// Diagonalization of :math:`H_{loc}`.
    atom_diag const &h_loc_diagonalization() const { return h_diag; }

    /// Histograms related to the performance analysis.
    histo_map_t const &get_performance_analysis() const { return _performance_analysis; }

    /// Monte Carlo average sign.
    mc_weight_t average_sign() const { return _average_sign; }

    /// Status of the ``solve()`` on exit.
    int solve_status() const { return _solve_status; }

    CPP2PY_IGNORE
    static std::string hdf5_format() { return "CTHYB_SolverCore"; }

    // Function that writes the solver_core to hdf5 file
    friend void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s) {
      h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);
      write_hdf5_format(grp, s);
      h5_write_attribute(grp, "TRIQS_GIT_HASH", std::string(STRINGIZE(TRIQS_GIT_HASH)));
      h5_write_attribute(grp, "CTHYB_GIT_HASH", std::string(STRINGIZE(CTHYB_GIT_HASH)));
      h5_write(grp, "container_set", s.container_set());
      h5_write(grp, "constr_parameters", s.constr_parameters);
      h5_write(grp, "solve_parameters", s.solve_parameters);
      h5_write(grp, "G0_iw", s._G0_iw);
      h5_write(grp, "Delta_tau", s._Delta_tau);

      h5_write(grp, "h_diag", s.h_diag);
      h5_write(grp, "h_loc", s._h_loc);
      h5_write(grp, "density_matrix", s._density_matrix);
      h5_write(grp, "average_sign", s._average_sign);
      h5_write(grp, "solve_status", s._solve_status);
      h5_write(grp, "Delta_infty_vec", s.Delta_infty_vec);
    }

    // Function that read all containers to hdf5 file
    CPP2PY_IGNORE
    static solver_core h5_read_construct(h5::group h5group, std::string subgroup_name) {
      h5::group grp          = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);
      auto constr_parameters = h5::h5_read<constr_parameters_t>(grp, "constr_parameters");
      auto s                 = solver_core{constr_parameters};
      h5_read(grp, "container_set", s.container_set());
      h5_read(grp, "solve_parameters", s.solve_parameters);
      h5_read(grp, "G0_iw", s._G0_iw);
      h5_read(grp, "Delta_tau", s._Delta_tau);

      h5_try_read(grp, "h_diag", s.h_diag);
      h5_try_read(grp, "h_loc", s._h_loc);
      h5_try_read(grp, "density_matrix", s._density_matrix);
      h5_try_read(grp, "average_sign", s._average_sign);
      h5_try_read(grp, "solve_status", s._solve_status);
      h5_try_read(grp, "Delta_infty_vec", s.Delta_infty_vec);

      return s;
    }
  };
} // namespace triqs_cthyb
