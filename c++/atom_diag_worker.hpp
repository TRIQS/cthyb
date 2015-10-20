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

// Division of Hilbert Space into sub hilbert spaces, using either autopartitioning or quantum numbers.
struct atom_diag_worker {

 atom_diag_worker(atom_diag* hdiag, int n_min = 0, int n_max = 10000) : hdiag(hdiag), n_min(n_min), n_max(n_max) {}

 void autopartition();
 void partition_with_qn(std::vector<many_body_op_t> const& qn_vector);

 private:
 atom_diag* hdiag;
 int n_min, n_max;
 
 // Create matrix of an operator acting from one subspace to another (returns matrix + number of its nonzero elements)
 matrix<double> make_op_matrix(imperative_operator<hilbert_space> const& op, int from_sp, int to_sp) const;

 void complete();
 bool fock_state_filter(fock_state_t s);
};
}
