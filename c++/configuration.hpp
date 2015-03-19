/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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
#pragma once
#include "sorted_spaces.hpp"
#include "./util.hpp"
#include <triqs/hilbert_space/hilbert_space.hpp>
#include <triqs/utility/time_pt.hpp>
#include <map>

namespace cthyb {

using triqs::utility::time_pt;
using triqs::utility::time_segment;

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

 configuration(double beta) : beta_(beta), id(0) {}

 double beta() const { return beta_; }
 int size() const { return oplist.size(); }

 void insert(time_pt tau, op_desc op) {oplist.insert({tau,op});}
 void erase(time_pt const & t) { oplist.erase(t); }

 oplist_t::iterator begin() { return oplist.begin(); }
 oplist_t::iterator end() { return oplist.end();}
 oplist_t::const_iterator begin() const { return oplist.begin();}
 oplist_t::const_iterator end() const { return oplist.end(); }
 
 friend std::ostream& operator<<(std::ostream& out, configuration const& c) {
  for (auto const& op : c)
   out << "tau = " << op.first << " : " << (op.second.dagger ? "Cdag(" : "C(") << op.second.block_index << ","
       << op.second.inner_index << ")\n";
  return out;
 }

 friend void h5_write(triqs::h5::group conf, std::string conf_group_name, configuration const& c) {
  triqs::h5::group conf_group = conf.create_group(conf_group_name);
  for (auto const& op : c) {
   // create group for given tau
   auto tau_group_name = std::to_string(double(op.first));
   triqs::h5::group tau_group = conf_group.create_group(tau_group_name);
   // in tau subgroup, write operator info
   h5_write(tau_group, "block", op.second.block_index);
   h5_write(tau_group, "inner", op.second.inner_index);
   h5_write(tau_group, "dagger", op.second.dagger);
  }
 }

 void print_to_h5(){
  std::string filename = "configs.h5";
  triqs::h5::file hfile(filename.c_str(), exists(filename) ? H5F_ACC_RDWR : H5F_ACC_TRUNC);
  h5_write(hfile, "c_"+std::to_string(this->id), *this);
  hfile.close();
 }

 template <class Archive> void serialize(Archive& ar, const unsigned int version) {
  ar& boost::serialization::make_nvp("oplist", oplist) & boost::serialization::make_nvp("beta", beta_);
 }

 int id; // configuration id, for debug purposes
 private:
 double beta_;
 oplist_t oplist;
};
}

