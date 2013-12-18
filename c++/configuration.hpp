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
 };

 // a map associating an operator to an imaginary time
 using oplist_t = std::map<time_pt, op_desc, std::greater<time_pt>>;
 // using oplist_t=boost::container::flat_map<time_pt, op_desc, std::greater<time_pt>> ;
 oplist_t oplist;

 // construction and the basics stuff. value semantics, except = ?
 configuration(double beta) : beta_(beta) {}

 double beta() const { return beta_; }

 friend std::ostream& operator<<(std::ostream& out, configuration const& c) {
  for (auto const& op : c.oplist)
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

