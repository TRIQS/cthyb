/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, H. U.R. Strand, M. Ferrero and O. Parcollet
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

#include <climits>

#include "./config.hpp"
#include "./types.hpp"

namespace triqs_cthyb {

  using namespace triqs::operators;
  using indices_map_t = std::map<triqs::operators::indices_t, triqs::operators::indices_t>;

  // All the arguments of the solver_core constructor
  struct constr_parameters_t {

    /// Inverse temperature
    double beta;

    ///block structure of the gf
    gf_struct_t gf_struct;

    /// Number of Matsubara frequencies for gf<imfreq, matrix_valued>
    int n_iw = 1025;

    /// Number of tau points for gf<imtime, matrix_valued>
    int n_tau = 10001;

    /// Number of Legendre polynomials for gf<legendre, matrix_valued>
    int n_l = 50;

    /// are Delta_tau and hloc_0 provided as input? if so they don't need to be calculated from G0_iw
    bool Delta_interface = false;

    /// Write constr_parameters_t to hdf5
    friend void h5_write(h5::group h5group, std::string subgroup_name, constr_parameters_t const &sp);

    /// Read constr_parameters_t from hdf5
    friend void h5_read(h5::group h5group, std::string subgroup_name, constr_parameters_t &sp);
  };

  // All the arguments of the solve function
  struct solve_parameters_t {

    /// Interacting part of the atomic Hamiltonian
    /// type: Operator
    many_body_op_t h_int;

    /// Number of QMC cycles
    int n_cycles;

    /// Partition method
    /// type: str
    std::string partition_method = "autopartition";

    /// Quantum numbers
    /// type: list(Operator)
    /// default: []
    std::vector<many_body_op_t> quantum_numbers = std::vector<many_body_op_t>{};

    /// Restrict local Hilbert space to states with at least this number of particles
    /// default: 0
    int loc_n_min = 0;

    /// Restrict local Hilbert space to states with at most this number of particles
    /// default: INT_MAX
    int loc_n_max = INT_MAX;

    /// Length of a single QMC cycle
    /// default: 50
    int length_cycle = 50;

    /// Number of cycles for thermalization
    /// default: 5000
    int n_warmup_cycles = 5000;

    /// Seed for random number generator
    /// default: 34788 + 928374 * MPI.rank
    int random_seed = 34788 + 928374 * mpi::communicator().rank();

    /// Name of random number generator
    /// type: str
    std::string random_name = "";

    /// Maximum runtime in seconds, use -1 to set infinite
    /// default: -1 = infinite
    int max_time = -1;

    /// Verbosity level
    /// default: 3 on MPI rank 0, 0 otherwise.
    int verbosity = ((mpi::communicator().rank() == 0) ? 3 : 0); // silence the slave nodes

    /// Add shifting an operator as a move?
    bool move_shift = true;

    /// Add double insertions as a move?
    bool move_double = true;

    /// Calculate the full trace or use an estimate?
    bool use_trace_estimator = false;

    /// Measure G(tau)? :math:`G_{ij}(\tau)=G_{ji}^*(\tau)` is enforced for the resulting G(tau)
    bool measure_G_tau = true;

    /// Measure G_l (Legendre)? Note, no hermiticity in G_l is ensured
    bool measure_G_l = false;

    /// Measure O_tau by insertion
    std::optional<std::pair<many_body_op_t, many_body_op_t>> measure_O_tau = std::optional<std::pair<many_body_op_t, many_body_op_t>>{};

    /// Minumum of operator insertions in: O_tau by insertion measure
    int measure_O_tau_min_ins = 10;

    /// Measure G^4(tau,tau',tau'') with three fermionic times.
    bool measure_G2_tau = false;

    /// Measure G^4(inu,inu',inu'') with three fermionic frequencies.
    bool measure_G2_iw = false;

    /// Measure G^4(inu,inu',inu'') with three fermionic frequencies.
    bool measure_G2_iw_nfft = false;

    /// Measure G^4(iomega,inu,inu') within the particle-particle channel.
    bool measure_G2_iw_pp = false;

    /// Measure G^4(iomega,inu,inu') within the particle-particle channel.
    bool measure_G2_iw_pp_nfft = false;

    /// Measure G^4(iomega,inu,inu') within the particle-hole channel.
    bool measure_G2_iw_ph = false;

    /// Measure G^4(iomega,inu,inu') within the particle-hole channel.
    bool measure_G2_iw_ph_nfft = false;

    /// Measure G^2(iomega,l,l') within the particle-particle channel.
    bool measure_G2_iwll_pp = false;

    /// Measure G^2(iomega,l,l') within the particle-hole channel.
    bool measure_G2_iwll_ph = false;

    /// Order of block indices in the definition of G^2.
    block_order measure_G2_block_order = block_order::AABB;

    /// List of block index pairs of G^2 to measure.
    /// default: measure all blocks
    std::set<std::pair<std::string, std::string>> measure_G2_blocks = (std::set<std::pair<std::string, std::string>>{});

    /// Number of imaginary time slices for G^4 measurement.
    int measure_G2_n_tau = 10;

    /// Number of bosonic Matsubara frequencies for G^4 measurement.
    int measure_G2_n_bosonic = 30;

    /// Number of fermionic Matsubara frequencies for G^4 measurement.
    int measure_G2_n_fermionic = 30;

    /// Number of Legendre coefficients for G^4(iomega,l,l') measurement.
    int measure_G2_n_l = 20;

    /// NFFT buffer size for G^4(iomega,l,l') measurement.
    int measure_G2_iwll_nfft_buf_size = 100;

    /// NFFT buffer sizes for different blocks
    /// default: 100 for every block
    std::map<std::string, long> nfft_buf_sizes = (std::map<std::string, long>{});

    /// Measure perturbation order?
    bool measure_pert_order = false;

    /// Measure the reduced impurity density matrix?
    bool measure_density_matrix = false;

    /// Use the norm of the density matrix in the weight if true, otherwise use Trace
    bool use_norm_as_weight = false;

    /// Analyse performance of trace computation with histograms (developers only)?
    bool performance_analysis = false;

    /// Operator insertion/removal probabilities for different blocks
    /// type: dict(str:float)
    /// default: {}
    std::map<std::string, double> proposal_prob = (std::map<std::string, double>{});

    /// List of global moves (with their names).
    /// Each move is specified with an index substitution dictionary.
    /// type: dict(str : dict(indices : indices))
    /// default: {}
    std::map<std::string, indices_map_t> move_global = (std::map<std::string, indices_map_t>{});

    /// Overall probability of the global moves
    double move_global_prob = 0.05;

    /// Threshold below which imaginary components of Delta and h_loc are set to zero
    double imag_threshold = 1.e-13;

    /// The maximum size of the determinant matrix before a resize
    int det_init_size = 100;

    /// Max number of ops before the test of deviation of the det, M^-1 is performed.
    int det_n_operations_before_check = 100;

    /// Threshold for determinant precision warnings
    double det_precision_warning = 1.e-8;

    /// Threshold for determinant precision error
    double det_precision_error = 1.e-5;

    /// Bound for the determinant matrix being singular, abs(det) > singular_threshold. If <0, it is !isnormal(abs(det))
    double det_singular_threshold = -1;

    solve_parameters_t() {}

    solve_parameters_t(many_body_op_t h_int, int n_cycles) : h_int(h_int), n_cycles(n_cycles) {}

    /// Write solve_parameters_t to hdf5
    friend void h5_write(h5::group h5group, std::string subgroup_name, solve_parameters_t const &sp);

    /// Read solve_parameters_t from hdf5
    friend void h5_read(h5::group h5group, std::string subgroup_name, solve_parameters_t &sp);

    /// Threshold below which which off diagonal components of hloc are set to 0
    double off_diag_threshold = 0.0;

    /// Quadratic part of the local Hamiltonian. Must be provided if the Delta interface is used
    std::optional<many_body_op_t> h_loc0 = {};
  };

  /// A struct combining both constr_params_t and solve_params_t
  struct params_t : constr_parameters_t, solve_parameters_t {
    params_t(constr_parameters_t constr_parameters_, solve_parameters_t solve_parameters_)
      : constr_parameters_t(constr_parameters_), solve_parameters_t(solve_parameters_) {}
  };
}
