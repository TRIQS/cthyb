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
#ifndef TRIQS_CTQMC_KRYLOV_MOVE_INSERT_H
#define TRIQS_CTQMC_KRYLOV_MOVE_INSERT_H
#include "./qmc_data.hpp"
#include <triqs/mc_tools/random_generator.hpp>
#include <triqs/mc_tools/histograms.hpp>

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

 typedef std::complex<double> mc_weight_type;

 //Insertion of C, C^dagger operator
 class move_insert_c_cdag {

  qmc_data & data;
  configuration & config;
  mc_tools::random_generator & rng;
  
  int block_index;
  int block_size;
  
  bool record_histograms;
  mc_tools::histogram_binned histo_length_proposed, histo_length_accepted;
  double delta_tau;
  qmc_data::trace_t new_trace;
  
  std::vector<configuration::oplist_t::iterator> inserted_ops;
  
  public :
  //-----------------------------------------------

  move_insert_c_cdag(int block_index, int block_size, qmc_data & data, mc_tools::random_generator & rng, bool record_histograms):
   data(data), config(data.config), rng(rng),
   block_index(block_index),
   block_size(block_size),
   record_histograms(record_histograms),
   histo_length_proposed (0,config.beta()),
   histo_length_accepted (0,config.beta())
  {}

  //---------------------

  mc_weight_type attempt() {

#ifdef KRYLOV_DEBUG
   std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
   std::cerr << "* Attempt for move_insert_c_cdag (block " << block_index << ")" << std::endl;
   std::cerr << "* Configuration before:" << std::endl;
   std::cerr << config;
#endif

   // Pick up the value of alpha and choose the operators
   configuration::op_desc op1 {block_index, rng(block_size), true},
                          op2 {block_index, rng(block_size), false};

   // Choice of times for insertion. Find the time as double and them put them on the grid.
   //time_pt tau1 (rng(config.beta()),config.beta());
   //time_pt tau2 (rng(config.beta()),config.beta());
   time_pt tau1 = time_pt::random (rng, config.beta(),config.beta());
   time_pt tau2 = time_pt::random (rng, config.beta(),config.beta());

#ifdef KRYLOV_DEBUG
   std::cerr << "* Proposing to insert:" << std::endl;
   std::cerr << "Cdag(" << op1.block_index << "," << op1.inner_index << ")";
   std::cerr << " at " << tau1 << std::endl;
   std::cerr << "C(" << op2.block_index << "," << op2.inner_index << ")";
   std::cerr << " at " << tau2 << std::endl;
#endif
   
   // record the length of the proposed insertion
   delta_tau = double(tau2-tau1);
   if (record_histograms) { histo_length_proposed << delta_tau; }

   // Insert the operators op1 and op2 at time tau1, tau2
   // 1- In the very exceptional case where the insert has failed because an operator is already sitting here
   // (cf std::map doc for insert return), we reject the move.
   // 2- If ok, we store the iterator to the inserted operators for later removal in reject if necessary
   auto r1 = config.oplist.insert (std::make_pair(tau1, op1));
   if (!r1.second) return 0;
   inserted_ops.push_back(r1.first);
   auto r2 = config.oplist.insert (std::make_pair(tau2, op2));
   if (!r2.second) return 0;
   inserted_ops.push_back(r2.first);
   
   new_trace = data.atomic_corr();
   auto trace_ratio = new_trace/data.trace;
   
   auto & det = data.dets[block_index];
   std::size_t det_size = det.size();
   
   // Find the position for insertion in the determinant   
   // NB : the determinant stores the C in decreasing time order.   
   std::size_t num_c_dag, num_c;
   for(num_c_dag=0; num_c_dag < det_size; ++num_c_dag) { if(det.get_x(num_c_dag).first < tau1) break; }
   for(num_c=0; num_c < det_size; ++num_c) { if(det.get_y(num_c).first < tau2) break; }
   
   // Insert in the det. Returns the ratio of dets (Cf det_manip doc).
   auto det_ratio = det.try_insert(num_c_dag, num_c,
                                   std::make_pair(tau1,op1.inner_index),
                                   std::make_pair(tau2,op2.inner_index));

   // acceptance probability
   mc_weight_type p = trace_ratio * det_ratio;
   double t_ratio = std::pow(block_size* config.beta() / double(det.size()+1), 2);

#ifdef KRYLOV_DEBUG
   std::cerr << "Trace ratio: " << trace_ratio << '\t';
   std::cerr << "Det ratio: " << det_ratio << '\t';
   std::cerr << "Prefactor: " << t_ratio << '\t';
   std::cerr << "Weight: " << p*t_ratio << std::endl;
   std::cerr << "* Configuration after: " << std::endl;
   std::cerr << config;
#endif
   
   return p*t_ratio;
  }

  //----------------

  mc_weight_type accept() {

#ifdef KRYLOV_DEBUG
   std::cerr << "* The move is accepted" << std::endl;
#endif

   inserted_ops.clear();
   data.dets[block_index].complete_operation();
   if (record_histograms) { histo_length_accepted << delta_tau; }

   data.update_sign();
   data.trace = new_trace;

#ifdef KRYLOV_DEBUG
   std::size_t sizee = data.dets[block_index].size();
   std::cerr << "Det x: " << std::endl;
   for(std::size_t i = 0; i < sizee; ++i) std::cerr << i << " -> " << data.dets[block_index].get_x(i).first << std::endl;
   std::cerr << "Det y: " << std::endl;
   for(std::size_t i = 0; i < sizee; ++i) std::cerr << i << " -> " << data.dets[block_index].get_y(i).first << std::endl;

   std::cerr << "Sign correction: " << data.current_sign / data.old_sign << std::endl;
   std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   return data.current_sign / data.old_sign;
  }

  //----------------

  void reject() {

#ifdef KRYLOV_DEBUG
   std::cerr << "* The move is rejected" << std::endl;
   std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif
   
   for (auto & it: inserted_ops) config.oplist.erase(it);
   inserted_ops.clear();
  }

 };

}}}}
#endif
