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
#include "./atom_diag.hpp"

namespace cthyb {

/// The atomic partition function
double partition_function(atom_diag const& atom, double beta);

/// The atomic density matrix
block_matrix_t atomic_density_matrix(atom_diag const& atom, double beta);

/// The atomic green function, possibly with excluded states (default none)
block_gf<imtime> atomic_gf(atom_diag const& atom, double beta, std::map<std::string, indices_t> const& indices_list, int n_tau,
                           std::vector<std::pair<int, int>> const& excluded_states = {});


/// Trace (op * density_matrix)
h_scalar_t trace_rho_op(block_matrix_t const& density_matrix, many_body_op_t const& op, atom_diag const& atom);

/// Act with operator op on state st
full_hilbert_space_state_t act(many_body_op_t const& op, full_hilbert_space_state_t const& st, atom_diag const& atom);

/** 
 * The operator op is supposed to be a quantum number (if not -> exception)
 * @return the eigenvalue by block
 */
std::vector<std::vector<double>> quantum_number_eigenvalues(many_body_op_t const& op, atom_diag const& atom);

/** 
 * The operator op is supposed to be a quantum number (if not -> exception)
 * @return the eigenvalue by block
 */
std::vector<std::vector<double>> quantum_number_eigenvalues2(many_body_op_t const& op, atom_diag const& atom);


}
