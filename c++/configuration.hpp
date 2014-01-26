#pragma once
#include "state.hpp"
#include "sorted_spaces.hpp"
#include "hilbert_space.hpp"
#include <triqs/utility/time_pt.hpp>
#include <map>

namespace cthyb_krylov {

using triqs::utility::time_pt;

// The configuration of the Monte Carlo
struct configuration {

 struct op_desc {  // The description of the C operator
  int block_index; // the block index of the operator
  int inner_index; // the inner index inside the block
  bool dagger;     // is the operator a dagger
  long linear_index;
  bool is_identity() const { return linear_index == -1;}
  static op_desc identity() { return {0,0,0,-1};} // -1 used for Id only
 };

 // a map associating an operator to an imaginary time
 using oplist_t = std::map<time_pt, op_desc, std::greater<time_pt>>;
 // using oplist_t=boost::container::flat_map<time_pt, op_desc, std::greater<time_pt>> ;
 private:
 oplist_t oplist;
 public : 

 // construction and the basics stuff. value semantics, except = ?
 configuration(double beta) : beta_(beta) {
  // add the two Id operators at the boundaries beta and 0
  oplist.insert( {time_pt::make_zero(beta), op_desc::identity()}); 
  oplist.insert( {time_pt::make_beta(beta), op_desc::identity()}); 
 }

 double beta() const { return beta_; }

 int size() const { return oplist.size() - 2; } // rm the 2 boundary points 

 auto insert(time_pt tau, op_desc op) DECL_AND_RETURN(oplist.insert({tau,op}));
 auto insert(std::pair<time_pt,op_desc> p) DECL_AND_RETURN(oplist.insert(p));

 void erase(oplist_t::iterator it) { oplist.erase(it); }

 oplist_t::iterator begin() { auto r = oplist.begin(); ++r; return r;}
 oplist_t::iterator end() { auto r = oplist.end(); --r; return r;}
 
 oplist_t::const_iterator begin() const { auto r = oplist.begin(); ++r; return r;}
 oplist_t::const_iterator end() const { auto r = oplist.end(); --r; return r;}
 
 oplist_t::const_iterator lowest_time_operator() const { auto r = end(); --r; return r;}
 oplist_t::const_iterator highest_time_operator() const { return begin();}
 oplist_t::const_iterator boundary_beta() const { return oplist.begin();}
 oplist_t::const_iterator boundary_zero() const { return end();}

 // returns the operator at the left of t (higher time). It can return the boundary Id operator at beta.
 oplist_t::const_iterator operator_just_after(time_pt t) const {
  auto r = oplist.lower_bound(t);
  return --r;
 }

 // returns the operator at the right of t (lower time). It can return the boundary Id operator at 0.
 oplist_t::const_iterator operator_just_before(time_pt t) const {
  return oplist.lower_bound(t);
 }

 friend std::ostream& operator<<(std::ostream& out, configuration const& c) {
  for (auto const& op : c)
   out << "tau = " << op.first << " : " << (op.second.dagger ? "Cdag(" : "C(") << op.second.block_index << ","
       << op.second.inner_index << ")\n";
  return out;
 }

 template <class Archive> void serialize(Archive& ar, const unsigned int version) {
  ar& boost::serialization::make_nvp("oplist", oplist) & boost::serialization::make_nvp("beta", beta_);
 }

 private:
 double beta_;
};
}

