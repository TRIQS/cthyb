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
#ifndef TRIQS_CTQMC_KRYLOV_ATOMIC_CORRELATOR_WORKER_H
#define TRIQS_CTQMC_KRYLOV_ATOMIC_CORRELATOR_WORKER_H
#include "configuration.hpp"
#include "sorted_spaces.hpp"
#include "exp_h_worker.hpp"
namespace cthyb_krylov {

 // A safe default if is_zero_state() isn't defined for StateType.
 template<typename StateType> bool is_zero_state(StateType const& st) { return false; }
    
 /**
  * A worker that computes the trace using krylov method, for a given configuration.
  * Has to live longer than the configuration...
  */
 class atomic_correlators_worker {

  public :
   
  typedef std::complex<double> result_t;
      
  private :
     
  const configuration * config; // must exists longer than this object.

  // The sorted space 
  sorted_spaces sosp;

  // the state
  typedef state<sub_hilbert_space,false> state_t;
  
  exp_h_worker<imperative_operator<sub_hilbert_space, false>, state_t> exp_h;

  std::vector<result_t> partial_traces; // the contribution of one block at the boundary
  result_t full_trace;
  
  // The minimal size of a matrix to be treated with exp_h_matrix
  std::size_t small_matrix_size;
  
  public : 
      
  atomic_correlators_worker(configuration & c, sorted_spaces const & sosp_,
                            double gs_energy_convergence, std::size_t small_matrix_size) :
    config(&c),
    sosp(sosp_),
    exp_h(sosp.hamiltonian(), sosp, gs_energy_convergence, small_matrix_size),
    partial_traces(c.boundary_block_states_ids.size(),0),
    small_matrix_size(small_matrix_size)
  {
  }

  // recompute and return the full trace
  result_t operator()() {
    full_trace = 0;

    for (int i=0; i < partial_traces.size(); ++i) {
      partial_traces[i] = compute_for_one_boundary_state(i);
      full_trace += partial_traces[i];
    }
    /*std::cout  << " number of BS "<< partial_traces.size() << std::endl;
    int i=0;
    for (auto x : partial_traces) if (std::abs(x)> 1.e-10) std::cout << i++ << "  "<< std::abs(x) << std::endl;
    std::cout  << "-----------"<<std::endl;
    */
    return full_trace;
  }
  
  // return the full trace, but recompute only one partial contribution
  result_t operator()(size_t n) {
      full_trace -= partial_traces[n];
      partial_traces[n] = compute_for_one_boundary_state(n);
      full_trace += partial_traces[n];
      return full_trace;
  }
  
  result_t compute_for_one_boundary_state(size_t n) {
   auto _begin = config->oplist.rbegin();
   auto _end   = config->oplist.rend();

   result_t trace = 0;
   
   std::size_t nsp, id;
   std::tie(nsp,id) = config->boundary_block_states_ids[n];
   
   state_t const& psi0 = sosp.get_eigensystems()[nsp].eigenstates[id];
       
   // do the first exp
   double dtau = ( _begin == _end ? config->beta() : double(_begin->first)); 
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

   //std::cout  << "looping "<< std::endl ;
   int cc =0;
   for (auto it = _begin; it != _end; cc++) { // do nothing if no operator

      // apply operator 
      auto const & op = sosp.get_fundamental_operator (it->second.dagger, it->second.block_index, it->second.inner_index);
      psi = op(psi);

      //std::cout << "gs energy interm "<< sosp.get_eigensystems()[psi.get_hilbert().get_index()].eigenvalues[0] << std::endl;
      // psi is already zero, makes no sense to proceed
      if(is_zero_state(psi)) { return 0;}
      //if(is_zero_state(psi)) { if (cc !=0) std::cout << "Cancel after : "<< cc << std::endl ; return 0;}
    
      // apply exponential. 
      double tau1 = double(it->first);
      ++it; 
      dtau = (it == _end ? config->beta() : double(it->first)) - tau1;  assert(dtau >0);        
      exp_h.apply (psi, dtau);
      //psi = exp_h (psi, dtau);
    }

    return dotc(psi0,psi);
  }

  const sorted_spaces& get_sorted_spaces() const { return sosp; };
  
 };

 //   template<class Archive> void serialize(Archive & ar, const unsigned int version) {
 //    ar & boost::make_nvp("oplist", oplist) & boost::make_nvp("beta",beta_) &  boost::make_nvp("boundary_block_state_ids",atomic_corr);
 //   }  

}
#endif

