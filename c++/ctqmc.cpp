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
#include "ctqmc.hpp"
#include <triqs/utility/callbacks.hpp>


#include "move_insert.hpp"
#include "move_remove.hpp"
#include "move_change_boundary_state.hpp"
#include "measure_z.hpp"
#include "measure_g.hpp"

#ifdef KRYLOV_STATS
#include "measure_boundary_state.hpp"
#endif

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

parameters ctqmc::apply_defaults(parameters const& p)
{
    parameter_defaults pdef;
        
    pdef.optional
        ("krylov_max_dim", krylov_params::default_max_dim, " unsigned int ")
        ("krylov_min_beta_threshold", krylov_params::default_min_beta_threshold, " double ")
        ("krylov_bs_use_cutoff", false, " bool ")
        ("krylov_bs_prob_cutoff", 1e-8, " double ") // put to <= 0 to include all boundary states.
#ifdef KRYLOV_STATS
        ("krylov_stats_file", std::string("krylov.stats.dat"), " string ")
        ("krylov_bs_stats_file", std::string("krylov.boundary_states.dat"), " string ")
#endif
    ;
    
    parameters updated_p = p;
    updated_p.update(pdef);
    return updated_p;
}
    
void ctqmc::solve() {

   // communicate on the world = all nodes
   boost::mpi::communicator c;

   // Moves
   auto & delta_names = data.delta.domain().names();
   for(size_t n = 0; n < data.delta.domain().size(); ++n){
     qmc.add_move(move_insert_c_cdag(n, data, qmc.rng(), false), "Insert Delta_" + delta_names[n]);
     qmc.add_move(move_remove_c_cdag(n, data, qmc.rng()), "Remove Delta_" + delta_names[n]);
   }
   
   if(!params["krylov_bs_use_cutoff"]){
        qmc.add_move(move_change_boundary_state(data, qmc.rng()), "Change the boundary state");
   }
   
   // Measurements
   qmc.add_measure(measure_z(), "Z measure");
   for(size_t n = 0; n < G_tau.domain().size(); ++n){
     qmc.add_measure(measure_g(n, G_tau[n], data), "G measure (" + G_tau.domain().names()[n] + ")");
   }
#ifdef KRYLOV_STATS
   if(!params["krylov_bs_use_cutoff"]){
       qmc.add_measure(measure_boundary_state(data, params["krylov_bs_stats_file"]), "Boundary state statistics");
   }
#endif

   // run!! The empty configuration has sign = 1
   qmc.start(1.0, triqs::utility::clock_callback(params["max_time"]));
   qmc.collect_results(c);
}

}}}}
