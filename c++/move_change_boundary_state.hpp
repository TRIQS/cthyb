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
#ifndef TRIQS_CTQMC_KRYLOV_MOVE_CHANGE_BOUNDARY_STATE_H
#define TRIQS_CTQMC_KRYLOV_MOVE_CHANGE_BOUNDARY_STATE_H
#include <cmath>
#include "./qmc_data.hpp"
#include <triqs/mc_tools/random_generator.hpp>

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

 typedef std::complex<double> mc_weight_type;

 /************************
   Insertion of C, C^dagger operator
  ****************************/
 
 class move_change_boundary_state {

  qmc_data & data;
  mc_tools::random_generator & rng;
  
  std::size_t subspace_index;
  std::size_t old_state_id;
  qmc_data::trace_t new_trace;
      
  public :
  //-----------------------------------------------

  move_change_boundary_state(qmc_data & data, mc_tools::random_generator & rng) :
    data(data),
    rng(rng)
  {}

  //---------------------

  mc_weight_type attempt()
  {

#ifdef EXT_DEBUG
   std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
   std::cerr << "* Attempt for move_change_boundary_state" << std::endl;
   std::cerr << "* Configuration before:" << std::endl;
   std::cerr << data.config;
#endif
      
    // Select a subspace
    subspace_index = rng(data.sosp.n_subspaces());
        
#ifdef EXT_DEBUG
   std::cerr << "* Proposing to change the boundry state in subspace " << subspace_index << std::endl;
#endif

    old_state_id = data.config.boundary_block_states_ids[subspace_index][0];
    data.config.boundary_block_states_ids[subspace_index][0] = rng(data.sosp.subspace(subspace_index).dimension());
    if(old_state_id == data.config.boundary_block_states_ids[subspace_index][0]) return 0;
    
    // Calculate the new trace (an update is made only for one block)
    new_trace = data.atomic_corr(subspace_index);
    auto ratio = new_trace/data.trace;
    
#ifdef EXT_DEBUG
   std::cerr << "Trace ratio: " << ratio << '\t';
   std::cerr << "Weight: " << ratio << std::endl;
   std::cerr << "* Configuration after: " << std::endl;
   std::cerr << data.config;
#endif
    
    return ratio;
  }

  //----------------

  mc_weight_type accept() {

#ifdef EXT_DEBUG
    std::cerr << "* The move is accepted" << std::endl;
    std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif      

    data.trace = new_trace;
      
    return 1.0;
  }

  //----------------

  void reject()
  {
#ifdef EXT_DEBUG
   std::cerr << "* The move is rejected" << std::endl;
   std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif
      
    // Restore the previously selected representative state
    data.config.boundary_block_states_ids[subspace_index][0] = old_state_id;
  }

 };

}}}}
#endif
