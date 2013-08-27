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
namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

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
  typedef state<partial_hilbert_space,false> state_t;
  
  exp_h_worker<imperative_operator<partial_hilbert_space, false>, state_t, true> exp_h_krylov;
  exp_h_worker<imperative_operator<partial_hilbert_space, false>, state_t, false> exp_h_matrix;

  std::vector<result_t> partial_traces; // the contribution of one block at the boundary
  result_t full_trace;
  
  // The minimal size of a matrix to be treated with exp_h_matrix
  std::size_t small_matrix_size;
  
  public : 
      
  atomic_correlators_worker(configuration & c, sorted_spaces const & sosp_, krylov_params kp, std::size_t small_matrix_size) :
    config(&c),
    sosp(sosp_),
    exp_h_krylov(sosp.hamiltonian(), kp), exp_h_matrix(),
    partial_traces(sosp.n_subspaces(),0),
    small_matrix_size(small_matrix_size)
  {
  }

  // recompute and return the full trace
  result_t operator()() {
    full_trace = 0;
    for (int i=0; i <sosp.n_subspaces(); ++i) {
      partial_traces[i] = compute_for_one_boundary_block(i);
      full_trace += partial_traces[i];
    }
    return full_trace;
  }
  
  // return the full trace, but recompute only one partial contribution
  result_t operator()(size_t bl) {
      full_trace -= partial_traces[bl];
      partial_traces[bl] = compute_for_one_boundary_block(bl);
      full_trace += partial_traces[bl];
      return full_trace;
  }
  
  result_t compute_for_one_boundary_block (size_t bl) {
   auto _begin = config->oplist.rbegin();
   auto _end   = config->oplist.rend();

   result_t trace = 0;
   
   auto const& eigensystem = sosp.get_eigensystems()[bl];
   auto const& eigenstates = eigensystem.eigenstates;
   
   // DEBUG
   bool use_krylov_worker = eigenstates.size() > small_matrix_size;
   
   for(std::size_t psi0_id : config->boundary_block_states_ids[bl]) {
     state_t const& psi0 = eigenstates[psi0_id];
       
     // do the first exp
     double dtau = ( _begin == _end ? config->beta() : double(_begin->first)); 
     state_t psi = eigenstates.size() > small_matrix_size ?
                   exp_h_krylov (psi0, dtau) :
                   exp_h_matrix (psi0, dtau, eigensystem);
     
     for (auto it = _begin; it != _end;) { // do nothing if no operator

        // apply operator 
        auto const & op = sosp.get_fundamental_operator (it->second.dagger, it->second.block_index, it->second.inner_index);
        psi = op(psi);
    
        // psi is already zero, makes no sense to proceed
        if(is_zero_state(psi)) goto next_psi0;
    
        // apply exponential. 
        double tau1 = double(it->first);
        ++it; 
        dtau = (it == _end ? config->beta() : double(it->first)) - tau1;  assert(dtau >0);
        
        auto const& psi_space = psi.get_hilbert();
        psi = psi_space.dimension() > small_matrix_size ?
              exp_h_krylov (psi, dtau) :
              exp_h_matrix (psi, dtau, sosp.get_eigensystems()[psi_space.get_index()]);
     }
     
     trace += dotc(psi0,psi);
     next_psi0:;
   }

   return trace;
  }

  const sorted_spaces& get_sorted_spaces() const { return sosp; };
  
 };

 //   template<class Archive> void serialize(Archive & ar, const unsigned int version) {
 //    ar & boost::make_nvp("oplist", oplist) & boost::make_nvp("beta",beta_) &  boost::make_nvp("boundary_block_state_ids",atomic_corr);
 //   }  

}}}}
#endif

