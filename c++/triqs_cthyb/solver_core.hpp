
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
#include <triqs/stat/histograms.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/atom_diag/functions.hpp>
#include <triqs/utility/macros.hpp>
#include <optional>

#include "types.hpp"
#include "container_set.hpp"
#include "parameters.hpp"
#include "configuration.hpp"

namespace triqs_cthyb {

  /// Continuous-time hybridization-expansion quantum Monte Carlo solver.
  class solver_core : public container_set_t {

    double beta;            // inverse temperature
    atom_diag h_diag;       // diagonalization of the local problem
    gf_struct_t gf_struct;  // Block structure of the Green function
    many_body_op_t _h_loc;  // The local Hamiltonian = h_int + h0
    many_body_op_t _h_loc0; //noninteracting part of the local Hamiltonian
    int n_iw, n_tau, n_l;
    bool delta_interface;

    std::vector<matrix_t> _density_matrix;            // density matrix, when used in Norm mode
    mpi::communicator _comm;                          // define the communicator, here MPI_COMM_WORLD
    histo_map_t _performance_analysis;                // Histograms used for performance analysis
    mc_weight_t _average_sign;                        // average sign of the QMC
    double _average_order;                            // average perturbation order
    double _auto_corr_time;                           // Auto-correlation time in units of MC cycles
    bool _auto_corr_time_converged = true;            // Whether the auto-correlation time estimate has saturated
    std::optional<std::map<std::string, nda::array<double, 1>>> _densities;        // Per-orbital densities
    std::optional<std::map<std::string, nda::array<double, 1>>> _densities_errors; // Density error bars
    int _solve_status;                                // Status of the solve upon exit: 0 for clean termination, > 0 otherwise.
    std::optional<configuration> _last_configuration; // Final configuration of the run

    // Single-particle Green's function containers
    std::optional<G_iw_t> _G0_iw;                                 // Non-interacting Matsubara Green's function
    G_tau_t _Delta_tau;                                           // Imaginary-time Hybridization function
    std::optional<std::vector<matrix<dcomplex>>> Delta_infty_vec; // Quadratic instantaneous part of G0_iw

    // Return reference to container_set
    container_set_t &container_set() { return static_cast<container_set_t &>(*this); }
    container_set_t const &container_set() const { return static_cast<container_set_t const &>(*this); }

    public:
    /// Parameters used for constructing the solver.
    constr_parameters_t constr_parameters;

    /// Parameters passed to the solve method.
    solve_parameters_t solve_parameters;

    /**
     * Construct a CTHYB solver.
     *
     * @param p Parameters used for constructing the solver.
     */
    solver_core(constr_parameters_t const &p);

    // Delete assignement operator because of const members
    solver_core(solver_core const &p)            = default;
    solver_core(solver_core &&p)                 = default;
    solver_core &operator=(solver_core const &p) = delete;
    solver_core &operator=(solver_core &&p)      = default;

    /**
     * Solve the impurity problem.
     *
     * @param p Parameters controlling the Monte Carlo simulation and measurements.
     */
    void solve(solve_parameters_t const &p);

    /// The local Hamiltonian \f$ H_{loc} \f$ used in the last solve.
    many_body_op_t h_loc() const { return _h_loc; }

    /// The noninteracting part of the local Hamiltonian.
    many_body_op_t h_loc0() const { return _h_loc0; }

    /// Parameters used for constructing the solver.
    constr_parameters_t last_constr_parameters() const { return constr_parameters; }

    /// Parameters used in the last solve.
    solve_parameters_t last_solve_parameters() const { return solve_parameters; }

    /// \f$ G_0^{-1}(i\omega_n = \infty) \f$ in Matsubara frequencies.
    [[deprecated("Use h_loc0() instead.")]]
    std::vector<matrix<dcomplex>> Delta_infty() {
      if (delta_interface) TRIQS_RUNTIME_ERROR << "Delta_infty cannot be accessed when using the Delta interface";
      return Delta_infty_vec.value();
    }

    /// Hybridization function \f$ \Delta(\tau) \f$ in imaginary time.
    block_gf_view<imtime> Delta_tau() { return _Delta_tau; }

    /// Non-interacting Green's function \f$ G_0(i\omega) \f$ in Matsubara frequencies.
    block_gf_view<imfreq> G0_iw() {
      if (delta_interface) TRIQS_RUNTIME_ERROR << "G0_iw cannot be accessed when using the Delta interface";
      return _G0_iw.value();
    }

    /// Atomic :math:`G(\tau)` in imaginary time.
    //block_gf_view<imtime> atomic_gf() const { return ::triqs_cthyb::atomic_gf(h_diag, beta, gf_struct, _Delta_tau[0].mesh().size()); }

    /// Accumulated density matrix.
    std::vector<matrix_t> density_matrix() const { return _density_matrix; }

    /// Diagonalization of \f$ H_{loc} \f$.
    atom_diag const &h_loc_diagonalization() const { return h_diag; }

    /// Histograms related to the performance analysis.
    C2PY_PROPERTY_GET(performance_analysis) histo_map_t get_performance_analysis() const { return _performance_analysis; }

    /// Monte Carlo average sign.
    mc_weight_t average_sign() const { return _average_sign; }

    /// Average perturbation order.
    double average_order() const { return _average_order; }

    /// Auto-correlation time in units of MC cycles.
    double auto_corr_time() const { return _auto_corr_time; }

    /// Whether the auto-correlation time estimate has saturated (false: it is only a lower bound, run longer).
    bool auto_corr_time_converged() const { return _auto_corr_time_converged; }

    /// Per-orbital densities, organized by blocks
    std::optional<std::map<std::string, nda::array<double, 1>>> densities() const { return _densities; }

    /// Error bars for per-orbital densities, organized by blocks
    std::optional<std::map<std::string, nda::array<double, 1>>> densities_errors() const { return _densities_errors; }

    /// Status of the ``solve()`` on exit.
    int solve_status() const { return _solve_status; }

    /// Final configuration of the last solve call.
    std::optional<configuration> last_configuration() const { return _last_configuration; }

    /// Is the solver compiled with support for complex hybridization?
    bool hybridisation_is_complex() const {
#ifdef HYBRIDISATION_IS_COMPLEX
      return true;
#else
      return false;
#endif
    }

    /// Is the solver compiled with support for a complex local Hamiltonian?
    bool local_hamiltonian_is_complex() const {
#ifdef LOCAL_HAMILTONIAN_IS_COMPLEX
      return true;
#else
      return false;
#endif
    }

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
      h5_write(grp, "average_order", s._average_order);
      h5_write(grp, "auto_corr_time", s._auto_corr_time);
      h5_write(grp, "auto_corr_time_converged", s._auto_corr_time_converged);
      h5_write(grp, "densities", s._densities);
      h5_write(grp, "densities_errors", s._densities_errors);
      h5_write(grp, "solve_status", s._solve_status);
      h5_write(grp, "Delta_infty_vec", s.Delta_infty_vec);
    }

    // Function that read all containers to hdf5 file
    C2PY_IGNORE static solver_core h5_read_construct(h5::group h5group, std::string subgroup_name) {
      h5::group grp          = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);
      auto constr_parameters = h5::h5_read<constr_parameters_t>(grp, "constr_parameters");
      auto s                 = solver_core{constr_parameters};
      h5_read(grp, "container_set", s.container_set());
      h5_read(grp, "solve_parameters", s.solve_parameters);
      h5_read(grp, "G0_iw", s._G0_iw);
      h5_read(grp, "Delta_tau", s._Delta_tau);

      h5::try_read(grp, "h_diag", s.h_diag);
      h5::try_read(grp, "h_loc", s._h_loc);
      h5::try_read(grp, "density_matrix", s._density_matrix);
      h5::try_read(grp, "average_sign", s._average_sign);
      h5::try_read(grp, "average_order", s._average_order);
      h5::try_read(grp, "auto_corr_time", s._auto_corr_time);
      h5::try_read(grp, "auto_corr_time_converged", s._auto_corr_time_converged);
      h5::try_read(grp, "densities", s._densities);
      h5::try_read(grp, "densities_errors", s._densities_errors);
      h5::try_read(grp, "solve_status", s._solve_status);
      h5::try_read(grp, "Delta_infty_vec", s.Delta_infty_vec);

      return s;
    }
  };
} // namespace triqs_cthyb
