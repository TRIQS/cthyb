#pragma once
#include <triqs/arrays.hpp>
//#include <ostream>
#include <unordered_map>
#include <cmath>
#include <boost/operators.hpp>
#include "hilbert_space.hpp"

namespace cthyb_krylov {

inline double conj(double x) { return x; } // questionable : collides with std::conj !!

// States of a Hilbert space : can either be described by a map
// or by a triqs::vector so there are two implementations controlled by BasedOnMap
template <typename HilbertSpace, typename ScalarType, bool BasedOnMap> class state {};

template <typename HilbertSpace, typename ScalarType, bool BasedOnMap>
state<HilbertSpace, ScalarType, BasedOnMap> make_zero_state(state<HilbertSpace, ScalarType, BasedOnMap> const& st) {
 return {st.get_hilbert()};
}

// -----------------------------------------------------------------------------------
// implementation based on a map : can work
// on huge hilbert spaces as long as there are not too
// many components in the state and not too many monomials
//  in the operator acting on the state...
// -----------------------------------------------------------------------------------
template <typename HilbertSpace, typename ScalarType>
class state<HilbertSpace, ScalarType, true> : boost::additive<state<HilbertSpace, ScalarType, true>>,
                                              boost::multiplicative<state<HilbertSpace, ScalarType, true>, ScalarType> {
 // derivations implement the vector space operations over ScalarType from the compounds operators +=, *=, ....
 const HilbertSpace* hs;
 using amplitude_t = std::unordered_map<std::size_t, ScalarType>;
 amplitude_t ampli;

 public:
 using value_type = ScalarType;

 state() : hs(nullptr) {} // non valid state !
 state(HilbertSpace const& hs_) : hs(&hs_) {}

 // Dimension of the Hilbert space
 // What if hs == nullptr ?
 int size() const { return hs->dimension(); }

 // Access to data
 value_type& operator()(int i) { return ampli[i]; }
 value_type const& operator()(int i) const { return ampli[i]; }

 // Basic operations
 state& operator+=(state const& s2) {
  for (auto const& aa : s2.ampli) {
   auto r = ampli.insert(aa);
   if (!r.second) r.first->second += aa.second;
  }
  prune();
  return *this;
 }

 state& operator-=(state const& s2) {
  for (auto const& aa : s2.ampli) {
   auto r = ampli.insert({aa.first, -aa.second});
   if (!r.second) r.first->second -= aa.second;
  }
  prune();
  return *this;
 }

 state& operator*=(value_type x) {
  for (auto& a : ampli) {
   a.second *= x;
  }
  prune();
  return *this;
 }

 state& operator/=(value_type x) { return operator*=(1 / x); }

 // Scalar product
 friend value_type dot_product(state const& s1, state const& s2) {
  value_type res = 0.0;
  for (auto const& a : s1.ampli) {
   if (s2.ampli.count(a.first) == 1) res += conj(a.second) * s2.ampli.at(a.first);
  }
  return res;
 }

 // Lambda (fs, amplitude)
 template<typename Lambda>
 friend void foreach(state const& st, Lambda l) {
  for (auto const& p : st.ampli) l(p.first, p.second);
 }

 //
 // Additions to StateVector concept
 //

 HilbertSpace const& get_hilbert() const { return *hs; }
 void set_hilbert(HilbertSpace const& hs_) { hs = &hs_; }

 private:
 void prune(double tolerance = 10e-10) {
  for (auto it = ampli.begin(); it != ampli.end(); it++) {
   if (std::fabs(it->second) < tolerance) ampli.erase(it);
  }
 }
};

// -----------------------------------------------------------------------------------
// implementation based on a vector
// -----------------------------------------------------------------------------------
template <typename HilbertSpace, typename ScalarType>
class state<HilbertSpace, ScalarType, false> : boost::additive<state<HilbertSpace, ScalarType, false>>,
                                               boost::multiplicative<state<HilbertSpace, ScalarType, false>, ScalarType> {

 const HilbertSpace* hs;
 using amplitude_t = triqs::arrays::vector<ScalarType>;
 amplitude_t ampli;

 public:
 using value_type = ScalarType;

 state() : hs(nullptr) {}
 state(HilbertSpace const& hs_) : hs(&hs_), ampli(hs_.dimension(), 0.0) {}

 // Dimension of the Hilbert space
 // What if hs == nullptr ?
 int size() const { return hs->dimension(); }

 // Access to data
 value_type& operator()(int i) { return ampli[i]; }
 value_type const& operator()(int i) const { return ampli[i]; }

 // Basic operations
 state& operator+=(state const& s2) {
  ampli += s2.ampli;
  return *this;
 }

 state& operator-=(state const& s2) {
  ampli -= s2.ampli;
  return *this;
 }

 state& operator*=(value_type x) {
  ampli *= x;
  return *this;
 }

 state& operator/=(value_type x) {
  ampli /= x;
  return *this;
 }

 // Scalar product
 friend value_type dot_product(state const& s1, state const& s2) { return dotc(s1.ampli, s2.ampli); }

 // Lambda (fs, amplitude)
 template<typename Lambda>
 friend void foreach(state const& st, Lambda l) {
  const auto L = st.size();
  for (size_t i = 0; i < L; ++i) l(i, st(i));
 }

 //
 // Additions to StateVector concept
 //

 // Full access to amplitudes
 amplitude_t const& amplitudes() const { return ampli; }
 amplitude_t& amplitudes() { return ampli; }

 // Get access to Hilbert space
 HilbertSpace const& get_hilbert() const { return *hs; }
 void set_hilbert(HilbertSpace const& hs_) { hs = &hs_; }

};

// Print state
template <typename HilbertSpace, typename ScalarType, bool BasedOnMap>
std::ostream& operator<<(std::ostream& os, state<HilbertSpace, ScalarType, BasedOnMap> const& s) {
 bool something_written = false;
 auto const& hs = s.get_hilbert();

 using value_type = typename state<HilbertSpace, ScalarType, BasedOnMap>::value_type;
 foreach(s, [&os,hs,&something_written](int i, value_type ampl){
  if (std::abs(ampl) > 1e-10){
   os << " +(" << ampl << ")" << "|" << hs.get_fock_state(i) << ">";
   something_written = true;
  }
 });

 if (!something_written) os << 0;
 return os;
}


}
