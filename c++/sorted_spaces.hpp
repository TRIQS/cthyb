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

#ifndef TRIQS_CTQMC_SORTED_SPACES
#define TRIQS_CTQMC_SORTED_SPACES

#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <memory>
#include <limits>
#include <algorithm> 
#include <triqs/utility/tuple_tools.hpp>
#include <triqs/arrays/vector.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/arrays/linalg/eigenelements.hpp>

#include "fundamental_operator_set.hpp"
#include "complete_hilbert_space.hpp"
#include "partial_hilbert_space.hpp"
#include "operator.hpp"
#include "imperative_operator.hpp"
#include "fock_state.hpp"
#include "state.hpp"

using namespace triqs::arrays;
using std::string;

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

// I define a more tolerant comparison between vectors
// for the quantum numbers
struct lt_dbl {
  bool operator()(std::vector<double> const & v1, std::vector<double> const & v2) const {
    for(int i=0; i<v1.size(); ++i) {
      if (v1[i] < (v2[i] - 1e-8)) return true;
      else if (v2[i] < (v1[i] - 1e-8)) return false;
    }
    return false;
  }
};

struct block_desc_t {
    std::string name;
    //std::vector<std::tuple<IndexType...>> indices;
    std::vector<fundamental_operator_set::indices_t> indices;
};

/*
  This class is used to divide the full Hilbert space into smaller
  subspaces using the quantum numbers.

*/
class sorted_spaces {

 using indices_t = fundamental_operator_set::indices_t;
 using one_indices_t = typename indices_t::value_type; // one element of the vector
  public:

  struct eigensystem_t {
      vector<double> eigenvalues;   // in ascending order
      std::vector<state<partial_hilbert_space,false>> eigenstates;
      matrix<double> unitary_matrix; // H = U * \Lambda * U^+
  };
      
  typedef double quantum_number_t;

  sorted_spaces(utility::many_body_operator<double> const& h_,
                std::vector<utility::many_body_operator<double>> const& qn_vector,
                fundamental_operator_set const& fops,
                std::vector<block_desc_t> const& block_structure):
    n_blocks(0), hamilt(h_, fops),
    creation_operators(fops.n_operators()), destruction_operators(fops.n_operators()),
    creation_map(fops.n_operators()), destruction_map(fops.n_operators()),
    creation_connection(fops.n_operators()), destruction_connection(fops.n_operators()) {

    std::map<indices_t, std::pair<int,int>> indices_to_ints;
    for(int bl=0; bl<block_structure.size(); ++bl){
        auto const& indices = block_structure[bl].indices;
        for(int i=0; i<indices.size(); ++i){
            indices_to_ints[indices[i]] = std::make_pair(bl,i);
        }
    }
        
    // create the map pair<int,int> --> int identifying operators
    for (auto const & ind_tuple: fops) {
      int_pair_to_n[ indices_to_ints.at(ind_tuple.first) ] = ind_tuple.second;
    }

    // the full Hilbert space
    complete_hilbert_space full_hs(fops);

    // make a vector of imperative operators for the quantum numbers
    for (auto & qn : qn_vector) {
      qn_operators.push_back( imperative_operator<complete_hilbert_space>(qn, fops) );
    }

    /*
      The first part consists in dividing the full Hilbert space
      into smaller subspaces using the quantum numbers
    */
    for (size_t r=0; r<full_hs.dimension(); ++r) {

      // fock_state corresponding to r
      fock_state fs = full_hs.get_fock_state(r);

      // the state we'll act on
      state<complete_hilbert_space, true> s(full_hs);
      s(r) = 1.0;

      // create the vector with the quantum numbers
      std::vector<quantum_number_t> qn = get_quantum_numbers(s);

      // if first time we meet these quantum numbers create partial Hilbert space
      if(map_qn_n.count(qn) == 0) {
        hilbert_spaces.push_back(std::make_shared<partial_hilbert_space>(hilbert_spaces.size()));
        quantum_numbers.push_back(qn);
        map_qn_n[qn] = n_blocks;
        n_blocks++;
      }

      // add fock state to partial Hilbert space
      hilbert_spaces[map_qn_n[qn]]->add_basis_fock(fs);
    }

    /*
      In this second part we want to derive the partial Hilbert space
      mapping. Basically we want to know if we act on a partial Hilbert
      space with a creation (destruction) operator, in which other
      partial Hilbert space we end up.
    */
    for (auto const & ind_tuple: fops) {

      // get the operators and their index
      int n = ind_tuple.second;
      auto create = utility::many_body_operator<double>::make_canonical(true, ind_tuple.first);
      auto destroy = utility::many_body_operator<double>::make_canonical(false, ind_tuple.first);
      //auto create = tuple::apply(triqs::utility::c_dag, ind_tuple.first);
      //auto destroy = tuple::apply(triqs::utility::c, ind_tuple.first);

      // construct their imperative counterpart
      imperative_operator<complete_hilbert_space> op_c_dag(create, fops);
      imperative_operator<complete_hilbert_space> op_c(destroy, fops);

      // to avoid declaring every time in the loop below
      std::vector<quantum_number_t> qn_before, qn_after;
      partial_hilbert_space * origin, * target;

      // these will be mapping tables
      creation_connection[n].resize(n_blocks, -1);
      destruction_connection[n].resize(n_blocks, -1);

      // now act on the state with the c, c_dag to see how quantum numbers change
      for (size_t r=0; r<full_hs.dimension(); ++r) {

        // the state we'll act on and its quantum numbers
        state<complete_hilbert_space, true> s(full_hs);
        s(r) = 1.0;
        qn_before = get_quantum_numbers(s);

        // apply creation on state to figure quantum numbers
        qn_after = get_quantum_numbers(op_c_dag(s));

        // insert in creation map checking whether it was already there
        if(dotc(op_c_dag(s), op_c_dag(s)) > 1.e-10) {
          origin = hilbert_spaces[map_qn_n[qn_before]].get();
          target = hilbert_spaces[map_qn_n[qn_after]].get();
          if (creation_map[n].count(origin) == 0) creation_map[n][origin] = target;
          else if (creation_map[n][origin] != target) std::cout << "error creation";
          creation_connection[n][map_qn_n[qn_before]] = map_qn_n[qn_after];
        }

        // apply destruction on state to figure quantum numbers
        qn_after = get_quantum_numbers(op_c(s));

        // insert in destruction map checking whether it was already there
        if(dotc(op_c(s), op_c(s)) > 1.e-10) {
          origin = hilbert_spaces[map_qn_n[qn_before]].get();
          target = hilbert_spaces[map_qn_n[qn_after]].get();
          if (destruction_map[n].count(origin) == 0) destruction_map[n][origin] = target;
          else if (destruction_map[n][origin] != target) std::cout << "error destruction";
          destruction_connection[n][map_qn_n[qn_before]] = map_qn_n[qn_after];
        }

      }

      // insert the creation and destruction operators in vectors. this is the fast version
      // of the operators because we explicitly use the map
      creation_operators[n] = imperative_operator<partial_hilbert_space, true>(create, fops, creation_map[n]);
      destruction_operators[n] = imperative_operator<partial_hilbert_space, true>(destroy, fops, destruction_map[n]);

    }
    
    // Compute energy levels and eigenvectors of the local Hamiltonian
    compute_eigensystems();
    
    // Shift the ground state energy of the local Hamiltonian to zero.
    for(auto & eigensystem : eigensystems) eigensystem.eigenvalues() -= get_gs_energy();
    hamilt = imperative_operator<partial_hilbert_space, false>(h_-get_gs_energy(), fops);
  }

  // get Hamiltonian
  imperative_operator<partial_hilbert_space, false> const & hamiltonian() const {
    return hamilt;
  }

  // get fundamental operators
  imperative_operator<partial_hilbert_space, true> const &
    get_fundamental_operator (bool dagger, int block_index, int inner_index) const {
      return (dagger ? creation_operators[int_pair_to_n.at(std::make_pair(block_index, inner_index))]
                     : destruction_operators[int_pair_to_n.at(std::make_pair(block_index, inner_index))]);
  }

  // connections for fundamental operators
  long fundamental_operator_connect(bool dagger, int block_index, int inner_index, size_t n) const {
    return (dagger ?  creation_connection[int_pair_to_n.at(std::make_pair(block_index, inner_index))][n]
                   :  destruction_connection[int_pair_to_n.at(std::make_pair(block_index, inner_index))][n]);
  }

  // subspaces
  partial_hilbert_space const & subspace(size_t n) const { return *hilbert_spaces[n]; }

  // substate
  state<partial_hilbert_space, false> substate(size_t n) const {
    return state<partial_hilbert_space, false>(*hilbert_spaces[n]);
  }

  // number of blocks
  size_t n_subspaces() const { return n_blocks; }

  // print some info
  friend std::ostream& operator<<(std::ostream& os, sorted_spaces const& ss) {

    os << "Number of blocks: " << ss.n_blocks << std::endl;
    for (int n=0; n<ss.hilbert_spaces.size(); ++n) {
      os << "Block " << n << ", ";
      os << "qn = ";
      std::for_each(ss.quantum_numbers[n].begin(), ss.quantum_numbers[n].end(), [&os](double x) { os << x << " "; });
      os << ", ";
      os << "size = " << ss.hilbert_spaces[n]->dimension() << std::endl;
    }
    return os;
  }
  
  std::vector<eigensystem_t> const& get_eigensystems() const { return eigensystems; }
  double get_gs_energy() const { return gs_energy; }

  private:

  typedef std::unordered_map<const partial_hilbert_space *, const partial_hilbert_space *> hilbert_map_t; 

  // Helper function to get quantum numbers
  std::vector<quantum_number_t> get_quantum_numbers(state<complete_hilbert_space> const & s) {
    std::vector<quantum_number_t> qn;
    for (auto & op : qn_operators) qn.push_back(dotc(s,op(s)));
    return qn;
  }
  
  void compute_eigensystems()
  {
    eigensystems.resize(n_subspaces());
    gs_energy = std::numeric_limits<double>::infinity();
    
    for(std::size_t spn=0; spn<n_subspaces(); ++spn){
        auto const& sp = subspace(spn);
        auto & eigensystem = eigensystems[spn];
        
        state<partial_hilbert_space,false> i_state(sp);
        matrix<double> h_matrix(sp.dimension(),sp.dimension()); 
        
        for(std::size_t i=0; i<sp.dimension(); ++i){
            i_state.amplitudes()() = 0;
            i_state(i) = 1;
            
            auto f_state = hamilt(i_state);
            h_matrix(range(),i)  = f_state.amplitudes();
         }
         linalg::eigenelements_worker<matrix_view<double>,true> ew(h_matrix);

         ew.invoke();
         eigensystem.eigenvalues = ew.values();
         eigensystem.unitary_matrix = h_matrix.transpose();
         gs_energy = std::min(gs_energy,eigensystem.eigenvalues[0]);
         
         eigensystem.eigenstates.reserve(sp.dimension());
         for(std::size_t e=0; e<sp.dimension(); ++e){
            eigensystem.eigenstates.emplace_back(sp);
            eigensystem.eigenstates.back().amplitudes() = h_matrix(e,range());
         }
    }    
  }
  
  // number of subspaces
  int n_blocks;

  // a map from a pair of integers to a single integer identifying the operator
  std::map<std::pair<int,int>, int> int_pair_to_n;

  // imperative operators for:
  // - the quantum number operators
  // - the hamiltonian,
  // - the creation and destruction operators
  std::vector<imperative_operator<complete_hilbert_space>> qn_operators;
  imperative_operator<partial_hilbert_space, false> hamilt;
  std::vector<imperative_operator<partial_hilbert_space, true>> creation_operators;
  std::vector<imperative_operator<partial_hilbert_space, true>> destruction_operators;

  // hilbert spaces and quantum numbers
  std::map<std::vector<double>, size_t, lt_dbl> map_qn_n;
  
  // have to use shared_ptr, because partial_hilbert_space is non-copyable
  std::vector<std::shared_ptr<partial_hilbert_space>> hilbert_spaces;
  
  std::vector<std::vector<quantum_number_t>> quantum_numbers;

  // mapping vectors
  std::vector<hilbert_map_t> creation_map;
  std::vector<hilbert_map_t> destruction_map;
  std::vector<std::vector<long>> creation_connection;
  std::vector<std::vector<long>> destruction_connection;
  
  // Eigensystem in each subspace
  std::vector<eigensystem_t> eigensystems;
  // Energy of the ground state
  double gs_energy;

};

}}}}

#endif
