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
#ifndef TRIQS_CTQMC_KRYLOV_OPERATOR_HPP_dfs87h
#define TRIQS_CTQMC_KRYLOV_OPERATOR_HPP_dfs87h

#include <ostream>
#include <istream>
#include <tuple>
#include <vector>
#include <map>
#include <limits>
#include <cmath>
#include <iterator>
#include <type_traits>
#include <boost/operators.hpp>
#include <boost/type_traits/has_less.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

#include "triqs/utility/tuple_tools.hpp"
#include "triqs/utility/tuple_serialize.hpp"
#include <triqs/utility/dressed_iterator.hpp>

namespace triqs { namespace utility {

 template<typename T> struct _meta_change_constchar_to_stdstring { using type = T;};
 template<> struct _meta_change_constchar_to_stdstring<const char *> { using type = std::string;};

 template<typename scalar_t, typename... IndexTypes>
  class many_body_operator :
   boost::addable<many_body_operator<scalar_t,IndexTypes...>,
   boost::subtractable<many_body_operator<scalar_t, IndexTypes...>,
   boost::multipliable<many_body_operator<scalar_t, IndexTypes...>,
   boost::addable2<many_body_operator<scalar_t, IndexTypes...>,      scalar_t,
   boost::subtractable2<many_body_operator<scalar_t, IndexTypes...>, scalar_t,
   boost::multipliable2<many_body_operator<scalar_t, IndexTypes...>, scalar_t
   >>>>>> 
 {
  static constexpr int n_indices = sizeof...(IndexTypes); 

  static constexpr scalar_t small_coeff = 100*std::numeric_limits<scalar_t>::epsilon();
  
  public:

  // The canonical operator : a dagger and some indices
  struct canonical_ops_t { 
   bool dagger; 
   using tuple_t = std::tuple<IndexTypes ...>;
   //using tuple_t = std::tuple<typename _meta_change_constchar_to_stdstring<IndexTypes>::type ...>;
   static_assert(boost::has_less<tuple_t>::value, "All indices must be LessThanComparable.");
   tuple_t indices;
   // dagger < non dagger, and then indices
   friend bool operator < (canonical_ops_t const & a, canonical_ops_t const & b) { return (a.dagger!=b.dagger ? a.dagger > b.dagger : (a.dagger ? a.indices < b.indices : a.indices > b.indices));}
   friend bool operator > (canonical_ops_t const & a, canonical_ops_t const & b) { return b<a; }
   friend bool operator== (canonical_ops_t const & a, canonical_ops_t const & b) { return (a.dagger==b.dagger && a.indices == b.indices);}
   template<class Archive> void serialize(Archive & ar, const unsigned int version) { ar & dagger & indices; }    
  };

  // Monomial: an ordered set of creation/annihilation operators and comparison
  typedef std::vector<canonical_ops_t> monomial_t;

  friend bool operator<(monomial_t const& m1, monomial_t const& m2) {
   return m1.size() != m2.size() ? m1.size() < m2.size() :
    std::lexicographical_compare(m1.begin(),m1.end(),m2.begin(),m2.end());
  }

  // Map of all monomials with coefficients
  typedef std::map<monomial_t,scalar_t> monomials_map_t;

  monomials_map_t monomials;

  many_body_operator() = default;
  many_body_operator(many_body_operator const&) = default;
  many_body_operator(many_body_operator &&) = default;
  many_body_operator& operator=(many_body_operator const &) = default;
#ifndef TRIQS_WORKAROUND_INTEL_COMPILER_BUGS
  many_body_operator& operator=(many_body_operator && o) = default;
#endif

  template<typename T> 
   many_body_operator(many_body_operator<T, IndexTypes...> const & x) { *this = x; }

  template<typename T> 
   many_body_operator & operator = (many_body_operator<T, IndexTypes...> const & x) {
    monomials.clear();
    for (auto const & y : x.monomials) monomials.insert( std::make_pair(y.first, y.second));
   }


  // We use utility::dressed_iterator to dress iterators 
  // _cdress is a simple struct of refs to dress the iterators (Cf doc)
  struct _cdress { 
   scalar_t coef; monomial_t const & monomial;
   _cdress(typename monomials_map_t::const_iterator _it) : monomial(_it->first), coef(_it->second){} 
  };
  typedef utility::dressed_iterator<typename monomials_map_t::const_iterator,_cdress>  const_iterator;

  public:

  // Iterators (only const!)
  const_iterator begin()  const noexcept { return monomials.begin(); }
  const_iterator end()    const noexcept { return monomials.end(); }
  const_iterator cbegin() const noexcept { return monomials.cbegin(); }
  const_iterator cend()   const noexcept { return monomials.cend(); }

  // Zero operator?
  bool is_zero() const { return monomials.empty(); }
  
  // Algebraic operations involving scalar_t constants
  many_body_operator operator-() const {
   auto tmp =*this;
   for(auto & m : tmp.monomials) m.second = -m.second;
   return tmp;
  }

  many_body_operator & operator+=(scalar_t alpha) {
   bool is_new_monomial;
   typename monomials_map_t::iterator it;
   std::tie(it,is_new_monomial) = monomials.insert(std::make_pair(monomial_t(0),alpha));
   if(!is_new_monomial){
    it->second += alpha;
    erase_zero_monomial(monomials,it);
   }
   return *this;
  }

  many_body_operator & operator-=(scalar_t alpha) { return (*this) += (-alpha);}

  friend many_body_operator operator-(scalar_t alpha, many_body_operator const& op) { return -op + alpha; }

  many_body_operator & operator*= (scalar_t alpha) {
   if(std::abs(alpha) < small_coeff){ monomials.clear(); } 
   else { for(auto & m : monomials) m.second *= alpha; }
   return *this;
  }

  // Algebraic operations
  many_body_operator & operator+=(many_body_operator const& op) {
   bool is_new_monomial;
   typename monomials_map_t::iterator it;
   for(auto const& m : op.monomials){
    std::tie(it,is_new_monomial) = monomials.insert(m);
    if(!is_new_monomial){
     it->second += m.second;
     erase_zero_monomial(monomials,it);
    }
   }
   return *this;
  }

  many_body_operator & operator-=(many_body_operator const& op) {
   bool is_new_monomial;
   typename monomials_map_t::iterator it;
   for(auto const& m : op.monomials){
    std::tie(it,is_new_monomial) = monomials.insert(std::make_pair(m.first,-m.second));
    if(!is_new_monomial){
     it->second -= m.second;
     erase_zero_monomial(monomials,it);
    }
   }
   return *this;
  }

  many_body_operator & operator*=(many_body_operator const& op) {
   monomials_map_t tmp_map; // product will be stored here
   for(auto const& m : monomials)
    for(auto const& op_m : op.monomials){
     // prepare an unnormalized product
     monomial_t product_m;
     product_m.reserve(m.first.size() + op_m.first.size());
     std::copy(m.first.begin(), m.first.end(), std::back_inserter(product_m));
     std::copy(op_m.first.begin(), op_m.first.end(), std::back_inserter(product_m));
     normalize_and_insert(product_m, m.second*op_m.second, tmp_map);
    }
   std::swap(monomials, tmp_map);
   return *this;
  }

  // Boost.Serialization
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) { ar & monomials; }    

  // Free factory functions
  template<typename... T> friend many_body_operator<double, typename _meta_change_constchar_to_stdstring<T>::type ...> c(T...);
  template<typename... T> friend many_body_operator<double, typename _meta_change_constchar_to_stdstring<T>::type ...> c_dag(T...);
  template<typename... T> friend many_body_operator<double, typename _meta_change_constchar_to_stdstring<T>::type ...> n(T...);

  private:

  // Normalize a monomial and insert into a map
  static void normalize_and_insert(monomial_t & m, scalar_t coeff, monomials_map_t & target)
  {
   // The normalization is done by employing a simple bubble sort algorithms.
   // Apart from sorting elements this function keeps track of the sign and
   // recursively calls itself if a permutation of two operators produces a new
   // monomial
   if(m.size() >= 2){
    bool is_swapped;
    do {
     is_swapped = false;
     for (std::size_t n = 1; n < m.size(); ++n){
      canonical_ops_t & prev_index = m[n-1];
      canonical_ops_t & cur_index = m[n];
      if(prev_index == cur_index) return;   // The monomial is effectively zero
      if(prev_index > cur_index){
       // Are we swapping C and C^+ with the same indices?
       // A bit ugly ...
       canonical_ops_t cur_index_flipped_type(cur_index);
       cur_index_flipped_type.dagger = !cur_index_flipped_type.dagger;
       if(prev_index == cur_index_flipped_type){
	monomial_t new_m;
	new_m.reserve(m.size() - 2);
	std::copy(m.begin(), m.begin() + n-1, std::back_inserter(new_m));
	std::copy(m.begin() + n+1, m.end(), std::back_inserter(new_m));

	normalize_and_insert(new_m, coeff, target);
       }
       coeff = -coeff;
       std::swap(prev_index, cur_index);
       is_swapped = true;
      }
     }
    } while(is_swapped);
   }

   // Insert the result
   bool is_new_monomial;
   typename monomials_map_t::iterator it;
   std::tie(it,is_new_monomial) = target.insert(std::make_pair(m, coeff));
   if(!is_new_monomial){
    it->second += coeff;
    erase_zero_monomial(target,it);
   }
  }

  // Erase a monomial with a close-to-zero coefficient.
  static void erase_zero_monomial(monomials_map_t & m, typename monomials_map_t::iterator & it) {
   if(std::abs(it->second) < small_coeff) m.erase(it);
  }

  friend std::ostream& operator<<(std::ostream &os, canonical_ops_t const & op) {
   if(op.dagger) os << "^+";
   return os << op.indices;
  }

  friend std::ostream& operator<<(std::ostream &os, monomial_t const& m) {
   for(auto const& c : m){ os << "C" << c; }
   return os;
  }

  // Print many_body_operator itself
  friend std::ostream& operator<<(std::ostream& os, many_body_operator const& op) {
   if(op.monomials.size() != 0){
    bool print_plus = false;
    for(auto const& m : op.monomials){
     os << (print_plus ? " + " : "" ) << m.second;
     if(m.first.size()) os << "*";
     os << m.first;
     print_plus = true;
    }
   } else
    os << "0";
   return os;
  }

 };

 // Free functions to make creation/annihilation operators
 template<typename... IndexTypes>
  many_body_operator<double, typename _meta_change_constchar_to_stdstring<IndexTypes>::type ...> c(IndexTypes... indices) {
  //many_body_operator<double, IndexTypes...> c(IndexTypes... indices) {
   typedef many_body_operator<double,typename _meta_change_constchar_to_stdstring<IndexTypes>::type ...> c_t;
   typedef typename c_t::canonical_ops_t canonical_ops_t;

   c_t tmp;
#ifndef TRIQS_WORKAROUND_INTEL_COMPILER_BUGS
   typename c_t::monomial_t m = {canonical_ops_t{false, std::make_tuple(indices...)}};
#else
   typename c_t::monomial_t m; m.push_back (canonical_ops_t{false, std::make_tuple(indices...)});
#endif
   tmp.monomials.insert(std::make_pair(m,1.0));
   return tmp;
  }

 template<typename... IndexTypes>
  many_body_operator<double,typename _meta_change_constchar_to_stdstring<IndexTypes>::type ...> c_dag(IndexTypes... indices) {
  //many_body_operator<double,IndexTypes...> c_dag(IndexTypes... indices) {
   typedef many_body_operator<double, typename _meta_change_constchar_to_stdstring<IndexTypes>::type ...> c_dag_t;
   typedef typename c_dag_t::canonical_ops_t canonical_ops_t;

   c_dag_t tmp;
#ifndef TRIQS_WORKAROUND_INTEL_COMPILER_BUGS
   typename c_dag_t::monomial_t m = {canonical_ops_t{true, std::make_tuple(indices...)}};
#else
   typename c_dag_t::monomial_t m; m.push_back (canonical_ops_t{true, std::make_tuple(indices...)});
#endif
   tmp.monomials.insert(std::make_pair(m,1.0));
   return tmp;    
  }

 template<typename... IndexTypes>
  many_body_operator<double,typename _meta_change_constchar_to_stdstring<IndexTypes>::type ...> n(IndexTypes... indices) {
  //many_body_operator<double,IndexTypes...> n(IndexTypes... indices) {
   typedef many_body_operator<double,typename _meta_change_constchar_to_stdstring<IndexTypes>::type ...> n_t;
   typedef typename n_t::canonical_ops_t canonical_ops_t;
   n_t tmp;
#ifndef TRIQS_WORKAROUND_INTEL_COMPILER_BUGS
   typename n_t::monomial_t m =
   {
    canonical_ops_t{true, std::make_tuple(indices...)},
    canonical_ops_t{false, std::make_tuple(indices...)}
   };
#else
   typename n_t::monomial_t m;
   m.push_back (canonical_ops_t{true, std::make_tuple(indices...)});
   m.push_back (canonical_ops_t{false, std::make_tuple(indices...)});
#endif
   tmp.monomials.insert(std::make_pair(m,1.0));

   return tmp;
  }

}}

#endif
