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

#ifndef TRIQS_CTQMC_KRYLOV_IMPERATIVE_OPERATOR
#define TRIQS_CTQMC_KRYLOV_IMPERATIVE_OPERATOR

#include <triqs/utility/exceptions.hpp>
#include <triqs/utility/tuple_tools.hpp>
#include "operator.hpp"
#include "fock_state.hpp"
#include "fundamental_operator_set.hpp"
#include <vector>

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

 // Reminder
 // fock state convention:
 // C^dag_0  ...  C^dag_k   | 0 >
 // operator convention:
 // C^dag_0 .. C^dag_i ... C^0  .. C^j

 // NOT HERE in fock_sate
//#define QUICK_APPLY_OP

 /*
    This class is the imperative version of the many_body_operator
    It is a class template on the Hilbert space type and has an optimization
    option UseMap which allows the user to give a map describing the connection
    between Hilbert space when acting on a state. This map is a map of *pointers*
    to the Hilbert spaces, to avoid unnessacry copies.

    If UseMap is false, the constructor takes two arguments:

    imperative_operator(many_body_op, fundamental_ops)

    a trivial identity connection map is then used.

    If UseMap is true, the constructor takes three arguments:

    imperative_operator(many_body_op, fundamental_ops, hilbert_map)
    */
 template<typename HilbertType, bool UseMap = false>
  class imperative_operator {

   // Pass it as a template parameter?
   typedef double scalar_t;

   // coeffs has the coefficient in front of the monomial
   // nondaggers has the destruction operators (left to right)
   // daggers has the construction operators (left to right)
   struct one_term { 
    scalar_t coeff;
#ifndef QUICK_APPLY_OP
    std::vector<int> daggers, nondaggers; 
#endif
    uint64_t d_mask, dag_mask, d_count_mask, dag_count_mask;
   };
   std::vector<one_term> all_terms;

   typedef std::unordered_map<const HilbertType *, const HilbertType *> hilbert_map_t;
   hilbert_map_t hilbert_map;

   public:

   size_t n_monomials() const {return all_terms.size();}

   imperative_operator() {}

   // constructor from a many_body_operator, a fundamental_operator_set and a map (UseMap = true)
   template<typename ...IndexTypes>
    imperative_operator(
      utility::many_body_operator<double,IndexTypes...> const & op,
      fundamental_operator_set<IndexTypes...> const & fops,
      hilbert_map_t hmap = hilbert_map_t() ){ 

     //std::cout  << " operator "<< op << std::endl ;
     hilbert_map = hmap;
     if ( (hilbert_map.size()==0) != !UseMap) TRIQS_RUNTIME_ERROR << "Internal error";

     // Reminder:
     // it->first is a vector< tuple<op_type, IndexTypes...> > == C^dag ... C^dag C ... C
     // it->second is the coefficient

     // The goal here is to have a transcription of the many_body_operator in terms
     // of simple vectors (maybe the code below could be more elegant)
     for (auto const & term : op) { 
      std::vector<int> dag, ndag;
      uint64_t d_mask=0, dag_mask=0;
      for (auto const & canonical_op : term.monomial ) { 
       (canonical_op.dagger ? dag : ndag).push_back(fops.index_tuple_to_linear(canonical_op.indices));
       (canonical_op.dagger ? dag_mask : d_mask) |= ( uint64_t(1) <<  fops.index_tuple_to_linear(canonical_op.indices));
      }
      auto compute_count_mask = [](std::vector<int> const & d) {
       //std::cout  << " cpmut mask "<< std::endl ;
       //for (auto x : d) std::cout  << x << std::endl ; 
       //auto dd = d;
       //std::sort(dd.begin(), dd.end());
       // I am using this assumption. If not general, I have to sort and keep the sign
       //if (d != dd) TRIQS_RUNTIME_ERROR << " Internal error : ndag dag not sorted ";
       uint64_t mask = 0;
       bool is_on = (d.size() %2 == 1);
       for (int i =0; i< 64; ++i) {
	if ( std::find(begin(d), end(d), i) != end(d) ) is_on = !is_on;
	else if (is_on) mask |= ( uint64_t(1) << i);
       }
       //std::cout  << "mask is "<< mask <<std::endl ;
       //std::cout  << "----------------------------------- "<< std::endl ;
       return mask;
      };
      uint64_t d_count_mask=compute_count_mask(ndag), dag_count_mask =compute_count_mask(dag);
      
      //if (!UseMap) std::cout  << all_terms.size() << " d_mask "<< d_mask << " size " << ndag.size() << "d_count_mask"<< d_count_mask <<std::endl ;
#ifdef QUICK_APPLY_OP
      all_terms.push_back( one_term {term.coef, d_mask, dag_mask, d_count_mask, dag_count_mask});
#else
      all_terms.push_back( one_term {term.coef, dag, ndag, d_mask, dag_mask, d_count_mask, dag_count_mask});
#endif
     }
    }

   // Regular type
   imperative_operator(imperative_operator const &) = default;
   imperative_operator & operator = (imperative_operator const &) = default;

   private : 
   template<typename StateType> 
    StateType get_target_st(StateType const & st, std::true_type use_map ) const {
     auto it = hilbert_map.find(&(st.get_hilbert()));
     if (it == hilbert_map.end()) return StateType();
     return StateType (*(it->second));
    }

   template<typename StateType>  
    StateType get_target_st(StateType const & st, std::false_type use_map ) const { return StateType (st.get_hilbert()); }

   bool parity_number_of_bits (uint64_t v) const { 
    // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetNaive    
/*    static const bool ParityTable256[256] = 
    {
#   define P2(n) n, n^1, n^1, n
#   define P4(n) P2(n), P2(n^1), P2(n^1), P2(n)
#   define P6(n) P4(n), P4(n^1), P4(n^1), P4(n)
     P6(0), P6(1), P6(1), P6(0)
    };
*/
    // OR, for 32-bit words:
    //v ^= v >> 16;
  
    // only ok until 16 orbitals ! assert this or put the >> 16  
    v ^= v >> 8; v ^= v >> 4; v ^= v >> 2; v ^= v >> 1; return v & 0x01;

    //v ^= v >> 8; return ParityTable256[v & 0xff];

    // 32 bits ONLY without this instruction 
    //v ^= v >> 32;
    // Variation:
    //unsigned char * p = (unsigned char *) &v;
    //return ParityTable256[p[0] ^ p[1] ^ p[2] ^ p[3]];

   }

   public:

   // act on a state and return a new state
   template<typename StateType>    
    StateType operator() (StateType const & st) const { 
     StateType target_st = get_target_st(st, std::integral_constant<bool, UseMap>()); //moved
     apply(st,target_st);
     return target_st;
    }
      
#ifdef QUICK_APPLY_OP
 
   template<typename StateType>    
    void apply (StateType const & st, StateType & target_st) const { 

     // loop over monomials
     for (int i=0; i<n_monomials(); ++i) {

      auto d_mask = all_terms[i].d_mask, dag_mask = all_terms[i].dag_mask;
      auto d_count_mask = all_terms[i].d_count_mask, dag_count_mask = all_terms[i].dag_count_mask;
      auto coef =  all_terms[i].coeff;

      foreach(st, [&]( uint64_t f2, typename StateType::data_t amplitude) { 

      //for (auto const & a: st) {

       //uint64_t f2 = a.fs;

       if ( (f2 & d_mask) != d_mask) return;
       f2 &= ~d_mask;

       if ( ( (f2 ^ dag_mask) & dag_mask ) != dag_mask) return;
       uint64_t f3 = ~( ~f2 & ~dag_mask );

       auto sign_is_minus = parity_number_of_bits((f2 & d_count_mask) ^ (f3 & dag_count_mask) );

       // update state vector in target Hilbert space
       auto ind = target_st.get_hilbert().get_state_index(f3);
       target_st(ind) += amplitude* coef * (sign_is_minus ? -1.0 : 1.0);
      });

     }
    }

#else


   template<typename StateType>    
    void apply (StateType const & st, StateType & target_st) const { 

     // loop over monomials
     for (int i=0; i<n_monomials(); ++i) {

      auto d_mask = all_terms[i].d_mask, dag_mask = all_terms[i].dag_mask;
      auto d_count_mask = all_terms[i].d_count_mask, dag_count_mask = all_terms[i].dag_count_mask;

      //uint64_t sign_acc =0;
      for (auto const & a: st) {

       int sign = 0;
       bool skip = false;


       fock_state newf = a.fs;

       try { 
	//#define DEBUG_TEST
#ifdef DEBUG_TEST

	uint64_t f2 = a.fs;
	auto s = 0, s1=0;
	bool skip1 = false, skip2 = false;

	if (d_mask !=0) { 
	 if ( (f2 & d_mask) != d_mask ) skip1=true;
	 else { 
	  f2 &= ~d_mask;
	  s = parity_number_of_bits( f2 & d_count_mask);
	  s1=s;
	 }}

	if (!skip1  &&  (dag_mask !=0)) { 
	 if ( ( (f2 ^ dag_mask) & dag_mask ) != dag_mask) skip2= true;
	 else { 
	  f2 = ~( ~f2 & ~dag_mask );
	  s += parity_number_of_bits( f2 & dag_count_mask);
	 }}

#endif

	// start from the right end of destruction operators
	for (auto n = all_terms[i].nondaggers.rbegin(); n != all_terms[i].nondaggers.rend(); ++n) {
	 if (newf.filling_number(*n) == 1) {
	  newf.set_to_0(*n);
	  sign += newf.num_particles_below(*n);
	 } else {
	  skip = true;
	  break;
	 }
	}

#ifdef DEBUG_TEST
	if (skip != skip1) TRIQS_RUNTIME_ERROR << " error skip1 " << a.fs << d_mask << dag_mask <<" "<< d_count_mask << dag_count_mask <<" "<<skip << skip1 ;
	if ( !skip && (sign%2 != s1%2)) TRIQS_RUNTIME_ERROR << " error sign1 " << a.fs << d_mask << dag_mask <<" "<< d_count_mask << dag_count_mask<<" " << sign << s1 ;
#endif

	if (skip) continue;

	// same for creation operators
	for (auto n = all_terms[i].daggers.rbegin(); n != all_terms[i].daggers.rend(); ++n) {
	 if (newf.filling_number(*n) == 0) {
	  newf.set_to_1(*n);
	  sign += newf.num_particles_below(*n);
	 } else {
	  skip = true;
	  break;
	 }
	}
#ifdef DEBUG_TEST
	if (skip != skip2) TRIQS_RUNTIME_ERROR << " error skip2 ";
	if (!skip && (sign%2 != s%2)) TRIQS_RUNTIME_ERROR << " error sign2 " << a.fs << d_mask << dag_mask <<" "<< sign<< s <<s1;
#endif
	if (skip) continue;

#ifdef DEBUG_TEST
	if (f2 != newf)   TRIQS_RUNTIME_ERROR <<  " != " << a.fs << " "<< f2 << " "<< newf ;;
	if (s %2 != sign %2)   TRIQS_RUNTIME_ERROR <<  " != sign  " << s << " "<< sign << " "<< newf ;;
#endif
       }
       catch( std::exception const & e) { std::cout  << e.what()<< std::endl ; throw;}

       sign = 1-2*(sign%2);

       // update state vector in target Hilbert space
       auto ind = target_st.get_hilbert().get_state_index(newf);
       target_st(ind) += a.amplitude* all_terms[i].coeff * sign;
      }
     }


     //return target_st;
    }
#endif

  };

}}}}
#endif
