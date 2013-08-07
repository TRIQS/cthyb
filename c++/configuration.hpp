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
#ifndef TRIQS_CTQMC_KRYLOV_CONFIGURATION_H
#define TRIQS_CTQMC_KRYLOV_CONFIGURATION_H

#include <map>
#include "state.hpp"
#include "sorted_spaces.hpp"
#include "partial_hilbert_space.hpp"
#include <triqs/utility/time_pt.hpp>

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

 using triqs::utility::time_pt;

 /**
   The configuration of the Monte Carlo
   */
 struct configuration {

  struct op_desc {     ///< The description of the C operator
   int block_index;    ///< the block index of the operator
   int inner_index;    ///< the inner index inside the block
   bool dagger;        ///< is the operator a dagger
  };

  // a map associating an operator to an imaginary time
  typedef std::map<time_pt, op_desc, std::greater<time_pt>> oplist_t;
  oplist_t oplist;

  // The boundary states
  std::vector<std::vector<std::size_t>> boundary_block_states_ids;

  // construction and the basics stuff. value semantics, except = ?
  configuration(double beta_) :
   beta_(beta_)
  {}
  
  void reset_boundary_block_states(sorted_spaces const & sosp)
  {
      boundary_block_states_ids.resize(sosp.n_subspaces());
      for(std::size_t nsp = 0; nsp < sosp.n_subspaces(); ++nsp){
          // Should initialize all boundary states with something nonzero
          boundary_block_states_ids[nsp].push_back(0);
      }
  }
  
  void fill_boundary_block_states(sorted_spaces const & sosp, double prob_tolerance = 0)
  {
      // Atomic partition function
      double z = 0;
      for(auto const& es : sosp.get_eigensystems())
          for(auto e : es.eigenvalues){
              z += exp(-beta_*(e-sosp.get_gs_energy()));
          }
              
      // Fill boundary_block_states with states which have significant weights
      boundary_block_states_ids.resize(sosp.n_subspaces());
      for(std::size_t nsp = 0; nsp < sosp.n_subspaces(); ++nsp){
          auto const& eigensystem = sosp.get_eigensystems()[nsp];
          for(size_t n=0; n<eigensystem.eigenvalues.size(); ++n){
              if(prob_tolerance == 0){
                  boundary_block_states_ids[nsp].push_back(n);
                  continue;
              }
              double prob = exp(-beta_*(eigensystem.eigenvalues[n]-sosp.get_gs_energy()))/z;
              if(prob>prob_tolerance) boundary_block_states_ids[nsp].push_back(n);
          }
      }
  }

  double beta() const {return beta_;}

  friend std::ostream & operator << (std::ostream & out, configuration const & c) {
    out << "boundary_block_states_ids:" << std::endl;
    for (size_t bl=0; bl< c.boundary_block_states_ids.size(); ++bl){
        out << "block " << bl << ":" << std::endl;
        for(auto const& st : c.boundary_block_states_ids[bl]) out << st << " ";
        out << std::endl;
    }
    out << std::endl;
    for (auto const& op : c.oplist)  
        out << "tau = "<< op.first << " : " << (op.second.dagger ? "Cdag(" : "C(") << op.second.block_index << ","<< op.second.inner_index<<")\n";
    return out;
  }

  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
   ar & boost::serialization::make_nvp("oplist", oplist) & boost::serialization::make_nvp("beta",beta_) 
    & boost::serialization::make_nvp("boundary_block_states_ids",boundary_block_states_ids);
  }
  
  private:

  double beta_;
 };

}}}}
#endif

