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
#include "atomic_correlators_worker.hpp"
#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>

namespace cthyb_krylov {

atomic_correlators_worker::atomic_correlators_worker(configuration& c, sorted_spaces const& sosp_, double gs_energy_convergence,
                                                     int small_matrix_size)
   : config(&c),
     sosp(sosp_),
     exp_h(sosp.get_hamiltonian(), sosp, gs_energy_convergence, small_matrix_size),
     small_matrix_size(small_matrix_size) {}

//------------------------------------------------------------------------------

atomic_correlators_worker::result_t atomic_correlators_worker::operator()() {
 result_t full_trace = 0;

 for (int n = 0; n < int(config->boundary_block_states_ids.size()); ++n) {

  auto _begin = config->oplist.rbegin();
  auto _end = config->oplist.rend();

  int nsp, id;
  std::tie(nsp, id) = config->boundary_block_states_ids[n];

  state_t const& psi0 = sosp.get_eigensystems()[nsp].eigenstates[id];

  // do the first exp
  double dtau = (_begin == _end ? config->beta() : double(_begin->first));
  state_t psi = exp_h(psi0, dtau);

  /*
   * // NEED TO CHANGE pointer to number with -1 when nothing....
   // first check of structural cancellation
   std::vector<int> blo(sosp.n_subspaces());
   auto n_blocks = sosp.n_subspaces();
   for (int u = 0; u < n_blocks; ++u) blo[u]=u;

   for (auto it = _begin; it != _end;) { // do nothing if no operator
    auto const & op = sosp.get_fundamental_operator (it->second.dagger, it->second.block_index, it->second.inner_index);
    bool all_zero = true;
    for (int u = 0; u < n_blocks; ++u) {
     if (blo[u] != -1) {
      all_zero = false;
      blo[u] = op.get_hilbert_connection(blo[u]);
     }
    }
    if (all_zero) return 0;
   }
 */

  // std::cout  << "looping "<< std::endl ;
  int cc = 0;
  bool psi_is_zero = false;
  for (auto it = _begin; it != _end; cc++) { // do nothing if no operator

   // apply operator
   auto const& op = sosp.get_fundamental_operator(it->second.dagger, it->second.block_index, it->second.inner_index);
   psi = op(psi);

   // std::cout << "gs energy interm "<< sosp.get_eigensystems()[psi.get_hilbert().get_index()].eigenvalues[0] << std::endl;
   // psi is already zero, makes no sense to proceed
   if (is_zero_state(psi)) {
    psi_is_zero = true;
    break;
   }
   // if(is_zero_state(psi)) { if (cc !=0) std::cout << "Cancel after : "<< cc << std::endl ; return 0;}

   // apply exponential.
   double tau1 = double(it->first);
   ++it;
   dtau = (it == _end ? config->beta() : double(it->first)) - tau1;
   assert(dtau > 0);
   exp_h.apply(psi, dtau);
   // psi = exp_h (psi, dtau);
  }

  if (psi_is_zero) continue;

  auto partial_trace = dot_product(psi0, psi);
  full_trace += partial_trace;

  // std::cout  << " partial trace "<< n << " "<< partial_trace<< std::endl;
  /*std::cout  << " number of BS "<< partial_traces.size() << std::endl;
  int i=0;
  for (auto x : partial_traces) if (std::abs(x)> 1.e-10) std::cout << i++ << "  "<< std::abs(x) << std::endl;
  std::cout  << "-----------"<<std::endl;
  */
 }
 return full_trace;
}
}
