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

#ifndef TRIQS_CTQMC_KRYLOV_STATE
#define TRIQS_CTQMC_KRYLOV_STATE

#include <map>
#include <ostream>
#include <boost/operators.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <tuple>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <triqs/arrays/vector.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>

#include "fock_state.hpp"

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

// Helper
template<typename T>
typename std::enable_if<boost::is_complex<T>::value, T>::type _conj(T && x) { return conj(std::forward<T>(x)); }
template<typename T>
typename std::enable_if<!boost::is_complex<T>::value, T>::type _conj(T && x) { return std::forward<T>(x); }


/*
  The class describing states of a Hilbert space.
  These states can either be described by a map
  or by a triqs::vector so there are two implementations
  controlled by the bool template argument BasedOnMap
*/
template<typename HilbertSpace, bool BasedOnMap = true>
class state {};

/*
  This is the implementation of a state based on a
  a map so that in principle we can work
  on huge hilbert spaces as long as there are not too
  many components in the state and not too many monomials
  in the operator acting on the state...
*/
template<typename HilbertSpace>
class state<HilbertSpace, true> : boost::addable< state<HilbertSpace, true>,
                                  boost::subtractable< state<HilbertSpace, true>,
                                  boost::dividable2< state<HilbertSpace, true>, double ,
                                  boost::multipliable2< state<HilbertSpace, true>, double > > > > {

  public:

    typedef double data_t;
    typedef data_t value_type;
    typedef std::unordered_map<std::size_t, data_t> amplitude_t;

  private:

    const HilbertSpace * hs;
    amplitude_t ampli;

  public:

    // empty constructor
    // This constructor creates an object which is logically 'just zero', irrespectively of any Hilbert space.
    // However, technically this object is invalid and must not be used in expressions. For the sake of convenience, 
    // I add a generic free function which detects this invalid state: is_zero_state()
    state(): hs(nullptr) {}

    // constructor on a Hilbert space
    state(HilbertSpace const & hs_): hs(&hs_) {}

    //std::size_t n_amplitudes() const { return ampli.size(); }
    amplitude_t const& amplitudes() const { return ampli; }
    amplitude_t & amplitudes() { return ampli; }
    
    friend std::size_t get_space_dim(state const & st)
    {
      assert(st.hs!=nullptr);
      return st.hs->dimension();
    }
    
    friend state make_zero_state(state const & st)
    {
        assert(st.hs!=nullptr);
        state zero_st(*st.hs);
        return zero_st;
    }

    data_t & operator()(std::size_t i) { return ampli[i]; }

    HilbertSpace const& get_hilbert() const { return *hs; }

    state(state const&) = default;
    state& operator=(state const&) = default;

    // iterator (concept satisfying)
    struct deref_struct { size_t index; fock_state fs; data_t amplitude; };
    class const_iterator: public boost::iterator_facade<const_iterator, deref_struct, boost::forward_traversal_tag, deref_struct> {
      public:
      const_iterator (amplitude_t::const_iterator it_, const state * st_): it(it_), st(st_) {}
      private:
      friend class boost::iterator_core_access;
      amplitude_t::const_iterator it;
      const state * st;
      void increment() { ++it; }
      //void decrement() { --it; }
      bool equal(const_iterator const & other) const { return(other.it == it); }
      deref_struct dereference() const { return {it->first, st->get_hilbert().get_fock_state(it->first), it->second}; }
    };
    const_iterator begin() const { return const_iterator(ampli.begin(), this); }
    const_iterator end() const { return const_iterator(ampli.end(), this); }

    // basic operations
    state& operator+=(state const& another_state) {
      bool new_amplitude;
      amplitude_t::iterator it;
      for(amplitude_t::const_reference aa : another_state.ampli){
        std::tie(it,new_amplitude) = ampli.insert(aa);
        if(!new_amplitude) it->second += aa.second;
      }
      prune();
      return *this;
    }

    state& operator-=(state const& another_state) {
      bool new_amplitude;
      amplitude_t::iterator it;
      for(amplitude_t::const_reference aa : another_state.ampli){
        std::tie(it,new_amplitude) = ampli.insert(std::make_pair(aa.first,-aa.second));
        if(!new_amplitude) it->second -= aa.second;
      }
      prune();
      return *this;
    }

    state& operator*=(data_t x) {
      for(auto &a : ampli) {
        a.second *= x;
      }
      prune();
      return *this;
    }

    state& operator/=(data_t x) {(*this) *= (1/x); return *this;}

    void prune(double tolerance = 10e-10) {
      for(auto it = ampli.begin(); it != ampli.end(); it++){
        if(std::fabs(it->second) < tolerance) ampli.erase(it);
      }
    }

    friend std::ostream& operator<<(std::ostream& os, state const& s) {
      for(amplitude_t::const_reference a : s.ampli){
        os << " +(" << a.second << ")" << s.hs->get_fock_state(a.first);
      }
      return os;
    }

    // scalar product
    friend data_t dotc(state const & s1, state const & s2) {
      data_t res = 0.0;
      for (amplitude_t::const_reference a: s1.ampli) {
        if (s2.ampli.count(a.first) == 1) res += _conj(a.second) * s2.ampli.at(a.first);
      }
      return res;
    }
    
    friend bool is_zero_state(state const& st)
    {
       return st.amplitudes().size() == 0;
    }

};

// Lambda (fs, amplitude)
template<typename HilbertSpace, typename Lambda>
void foreach (state<HilbertSpace, true> const & st, Lambda l) { 
 for ( auto const & p : st.amplitudes() ) l( st.get_hilbert().get_fock_state(p.first) , p.second);
} 




/*
   This is the implementation of a state based on a
   a triqs::vector.
   */
template<typename HilbertSpace>
class state<HilbertSpace, false> : boost::addable< state<HilbertSpace, false>,
      boost::subtractable< state<HilbertSpace, false>,
      boost::dividable2< state<HilbertSpace, false>, double,
      boost::multipliable2< state<HilbertSpace, false>, double > > > > {


       public:

	typedef double data_t;
	typedef data_t value_type;
	typedef triqs::arrays::vector<data_t> amplitude_t;

       private:

	const HilbertSpace * hs;
	amplitude_t ampli;

       public:

	// default constructor
	state(): hs(nullptr) {}

	// constructor from just a hilbert_space
	state(HilbertSpace const & hs_): hs(&hs_), ampli(hs_.dimension(), 0.0) {}

	// value
	state(state const &) = default;
	state(state &&) = default;
	state & operator = (state const &) = default;

	// how many amplitudes
	size_t n_amplitudes() const { return ampli.size(); }

	friend std::size_t get_space_dim(state const & st)
	{
	 assert(st.hs!=nullptr);
	 return st.hs->dimension();
	}

	friend state make_zero_state(state const & st)
    {
        assert(st.hs!=nullptr);
        state zero_st(*st.hs);
        zero_st.amplitudes()() = 0;
        return zero_st;
    }
	
	// full access to amplitudes
	amplitude_t const & amplitudes() const { return ampli; }
	amplitude_t & amplitudes() { return ampli; }

	// access to data
	data_t       & operator()(std::size_t i) { return ampli[i]; }
	data_t const & operator()(std::size_t i) const { return ampli[i]; }

	// get access to hilbert space
	HilbertSpace const & get_hilbert() const { return *hs; }

	// iterator (concept satisfying)
	struct deref_struct { int index; fock_state fs; data_t amplitude; };
	struct const_iterator: public boost::iterator_facade<const_iterator, deref_struct, boost::bidirectional_traversal_tag, deref_struct> {
	 public :
	  const_iterator (int p_, const state * st_): p(p_), st(st_) {}
	 private:
	  friend class boost::iterator_core_access;
	  int p;
	  const state * st;
	  void increment() { ++p; }
	  void decrement() { --p; }
	  bool equal(const_iterator const & other) const { return(other.p == p); }
	  deref_struct dereference() const { return {p, st->get_hilbert().get_fock_state(p), st->ampli[p]}; }
	};
	const_iterator begin() const { return const_iterator(0, this); }
	const_iterator end() const { return const_iterator(ampli.size(), this); }

	// print
	friend std::ostream& operator<<(std::ostream& os, state const& s) {
	 for(int i=0; i<s.n_amplitudes(); ++i) {
      auto ampl = s(i);
      if(std::abs(ampl)<1e-10) continue;
	  os << " +(" << ampl << ")" << s.hs->get_fock_state(i);
	 }
	 return os;
	}

	// basic operations
	state & operator += (state const & another_state) {
	 ampli += another_state.ampli;
	 return *this;
	}

	state & operator -= (state const & another_state) {
	 ampli -= another_state.ampli;
	 return *this;
	}

	state & operator *= (data_t x) {
	 ampli *= x;
	 return *this;
	}

	state & operator /= (data_t x) {
	 ampli /= x;
	 return *this;
	}

	// scalar product
	friend data_t dotc(state const & s1, state const & s2) {
	 return dotc(s1.ampli, s2.ampli);
	}

	friend bool is_zero_state(state const& st, double tolerance = 1e-18)
	{
	 if (st.amplitudes().size() == 0) return true;
	 for (auto const &a : st.amplitudes()) if(std::fabs(a) > tolerance) return false;
	 return true;
	}

      };


// Lambda (fs, amplitude)
template<typename HilbertSpace, typename Lambda>
void foreach (state<HilbertSpace, false> const & st, Lambda l) { 
 const auto L = st.amplitudes().size();
 for (size_t i =0; i<L; ++i) l(st.get_hilbert().get_fock_state(i), st.amplitudes()[i]);
} 


}}}}
#endif
