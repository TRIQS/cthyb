/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by I. Krivenko, M. Ferrero, O. Parcollet
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
#include "configuration.hpp"
#include "sorted_spaces.hpp"
#include "exp_h_worker.hpp"
namespace cthyb_krylov {

/**
 * A worker that computes the trace using krylov method, for a given configuration.
 * Has to live longer than the configuration...
 */
class atomic_correlators_worker {
 public:
 
 using result_t = std::complex<double>;
 
 atomic_correlators_worker(configuration& c, sorted_spaces const& sosp_, double gs_energy_convergence, int small_matrix_size);

 result_t operator()(); // recompute and return the full trace

 const sorted_spaces& get_sorted_spaces() const { return sosp; }
 
 private:
 const configuration* config; // must exists longer than this object.
 sorted_spaces sosp;          // The sorted space
 int small_matrix_size;// The minimal size of a matrix to be treated with exp_h_matrix

 using state_t = state<sub_hilbert_space, double, false>;
 exp_h_worker<imperative_operator<sub_hilbert_space, false>, state_t> exp_h;
};
}
