#pragma once
#include "state.hpp"
#include "sorted_spaces.hpp"
#include "hilbert_space.hpp"
#include <triqs/utility/time_pt.hpp>
#include <map>

namespace cthyb_matrix {

using triqs::utility::time_pt;

struct op_desc {  // The description of the C operator
 int block_index; // the block index of the operator
 int inner_index; // the inner index inside the block
 bool dagger;     // is the operator a dagger
 long linear_index;
};

// The configuration of the Monte Carlo
struct configuration {

 // a map associating an operator to an imaginary time
 using oplist_t = std::map<time_pt, op_desc, std::greater<time_pt>>;

 configuration(double beta) : beta_(beta) {}

 double beta() const { return beta_; }
 int size() const { return oplist.size(); }

 void insert(time_pt tau, op_desc op) {oplist.insert({tau,op});}
 void erase(time_pt const & t) { oplist.erase(t); }

 oplist_t::iterator begin() { return oplist.begin(); }
 oplist_t::iterator end() { return oplist.end();}
 oplist_t::const_iterator begin() const { return oplist.begin();}
 oplist_t::const_iterator end() const { return oplist.end(); }
 
 friend std::ostream& operator<<(std::ostream& out, configuration const& c) {
  out << " Config Id = "<< c.id<< std::endl;
  for (auto const& op : c)
   out << "tau = " << op.first << " : " << (op.second.dagger ? "Cdag(" : "C(") << op.second.block_index << ","
       << op.second.inner_index << ")\n";
  return out;
 }

 template <class Archive> void serialize(Archive& ar, const unsigned int version) {
  ar& boost::serialization::make_nvp("oplist", oplist) & boost::serialization::make_nvp("beta", beta_);
 }

 uint64_t id = 0;
 bool print_debug() const {
  return false;
  int n = 3743;
  //return true; 
   if (id > n) throw "";
  return (id >= n) && (id < n + 1);
 }
 private:
 double beta_;
 oplist_t oplist;
};
}

