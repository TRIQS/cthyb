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

#include "./parameters.hpp"

namespace cthyb {

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, constr_parameters_t const &cp) {
    triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);
    h5_write(grp, "beta", cp.beta);
    h5_write(grp, "gf_struct", cp.gf_struct);
    h5_write(grp, "n_iw", cp.n_iw);
    h5_write(grp, "n_tau", cp.n_tau);
    h5_write(grp, "n_l", cp.n_l);
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, constr_parameters_t &cp) {
    triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);
    h5_read(grp, "beta", cp.beta);
    h5_read(grp, "gf_struct", cp.gf_struct);
    h5_read(grp, "n_iw", cp.n_iw);
    h5_read(grp, "n_tau", cp.n_tau);
    h5_read(grp, "n_l", cp.n_l);
  }

  void h5_write(triqs::h5::group h5group, std::string subgroup_name, solve_parameters_t const &sp) {
    triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);
    h5_write(grp, "h_int", sp.h_int);

    h5_write(grp, "n_cycles", sp.n_cycles);
    h5_write(grp, "partition_method", sp.partition_method);
    h5_write(grp, "quantum_numbers", sp.quantum_numbers);

    h5_write(grp, "length_cycle", sp.length_cycle);
    h5_write(grp, "n_warmup_cycles", sp.n_warmup_cycles);
    h5_write(grp, "random_seed", sp.random_seed);
    h5_write(grp, "random_name", sp.random_name);
    h5_write(grp, "max_time", sp.max_time);
    h5_write(grp, "verbosity", sp.verbosity);

    h5_write(grp, "move_shift", sp.move_shift);
    h5_write(grp, "move_double", sp.move_double);
    h5_write(grp, "use_trace_estimator", sp.use_trace_estimator);

    h5_write(grp, "measure_G_tau", sp.measure_G_tau);
    h5_write(grp, "measure_G_l", sp.measure_G_l);
    h5_write(grp, "measure_G2_tau", sp.measure_G2_tau);
    h5_write(grp, "measure_G2_iw", sp.measure_G2_iw);
    h5_write(grp, "measure_G2_iw_pp", sp.measure_G2_iw_pp);
    h5_write(grp, "measure_G2_iw_ph", sp.measure_G2_iw_ph);
    h5_write(grp, "measure_G2_iwll_pp", sp.measure_G2_iwll_pp);
    h5_write(grp, "measure_G2_iwll_ph", sp.measure_G2_iwll_ph);
    h5_write(grp, "measure_G2_block_order", sp.measure_G2_block_order);
    h5_write(grp, "measure_G2_blocks", sp.measure_G2_blocks);

    h5_write(grp, "measure_G2_n_tau", sp.measure_G2_n_tau);
    h5_write(grp, "measure_G2_n_bosonic", sp.measure_G2_n_bosonic);
    h5_write(grp, "measure_G2_n_fermionic", sp.measure_G2_n_fermionic);
    h5_write(grp, "measure_G2_n_l", sp.measure_G2_n_l);

    h5_write(grp, "measure_G2_iwll_nfft_buf_size", sp.measure_G2_iwll_nfft_buf_size);
    h5_write(grp, "nfft_buf_sizes", sp.nfft_buf_sizes);

    h5_write(grp, "measure_pert_order", sp.measure_pert_order);
    h5_write(grp, "measure_density_matrix", sp.measure_density_matrix);
    h5_write(grp, "use_norm_as_weight", sp.use_norm_as_weight);
    h5_write(grp, "performance_analysis", sp.performance_analysis);
    h5_write(grp, "proposal_prob", sp.proposal_prob);

    h5_write(grp, "move_global", sp.move_global);
    h5_write(grp, "move_global_prob", sp.move_global_prob);

    h5_write(grp, "imag_threshold", sp.imag_threshold);
  }

  void h5_read(triqs::h5::group h5group, std::string subgroup_name, solve_parameters_t &sp) {
    triqs::h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);
    h5_read(grp, "h_int", sp.h_int);

    h5_read(grp, "n_cycles", sp.n_cycles);
    h5_read(grp, "partition_method", sp.partition_method);
    h5_read(grp, "quantum_numbers", sp.quantum_numbers);

    h5_read(grp, "length_cycle", sp.length_cycle);
    h5_read(grp, "n_warmup_cycles", sp.n_warmup_cycles);
    h5_read(grp, "random_seed", sp.random_seed);
    h5_read(grp, "random_name", sp.random_name);
    h5_read(grp, "max_time", sp.max_time);
    h5_read(grp, "verbosity", sp.verbosity);

    h5_read(grp, "move_shift", sp.move_shift);
    h5_read(grp, "move_double", sp.move_double);
    h5_read(grp, "use_trace_estimator", sp.use_trace_estimator);

    h5_read(grp, "measure_G_tau", sp.measure_G_tau);
    h5_read(grp, "measure_G_l", sp.measure_G_l);
    h5_read(grp, "measure_G2_tau", sp.measure_G2_tau);
    h5_read(grp, "measure_G2_iw", sp.measure_G2_iw);
    h5_read(grp, "measure_G2_iw_pp", sp.measure_G2_iw_pp);
    h5_read(grp, "measure_G2_iw_ph", sp.measure_G2_iw_ph);
    h5_read(grp, "measure_G2_iwll_pp", sp.measure_G2_iwll_pp);
    h5_read(grp, "measure_G2_iwll_ph", sp.measure_G2_iwll_ph);
    h5_read(grp, "measure_G2_block_order", sp.measure_G2_block_order);
    h5_read(grp, "measure_G2_blocks", sp.measure_G2_blocks);

    h5_read(grp, "measure_G2_n_tau", sp.measure_G2_n_tau);
    h5_read(grp, "measure_G2_n_bosonic", sp.measure_G2_n_bosonic);
    h5_read(grp, "measure_G2_n_fermionic", sp.measure_G2_n_fermionic);
    h5_read(grp, "measure_G2_n_l", sp.measure_G2_n_l);

    h5_read(grp, "measure_G2_iwll_nfft_buf_size", sp.measure_G2_iwll_nfft_buf_size);
    h5_read(grp, "nfft_buf_sizes", sp.nfft_buf_sizes);

    h5_read(grp, "measure_pert_order", sp.measure_pert_order);
    h5_read(grp, "measure_density_matrix", sp.measure_density_matrix);
    h5_read(grp, "use_norm_as_weight", sp.use_norm_as_weight);
    h5_read(grp, "performance_analysis", sp.performance_analysis);
    h5_read(grp, "proposal_prob", sp.proposal_prob);

    h5_read(grp, "move_global", sp.move_global);
    h5_read(grp, "move_global_prob", sp.move_global_prob);

    h5_read(grp, "imag_threshold", sp.imag_threshold);
  }
  
} // namespace cthyb
