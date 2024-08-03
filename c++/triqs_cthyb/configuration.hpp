/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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
#include "./util.hpp"
#include <triqs/hilbert_space/hilbert_space.hpp>
#include <triqs/utility/time_pt.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/atom_diag/functions.hpp>

#include <h5/h5.hpp>

#include <map>

namespace triqs_cthyb {

  using triqs::utility::time_pt;
  using triqs::utility::time_segment;

  // The description of the C operator
  struct op_desc {
    int block_index;   // the block index of the operator
    int inner_index;   // the inner index inside the block
    bool dagger;       // is the operator a dagger
    long linear_index; // the cumulative index

    friend std::ostream &operator<<(std::ostream &out, op_desc const &op) {
      out << (op.dagger ? "Cdag(" : "C(") << op.block_index << "," << op.inner_index << ")";
      return out;
    }

    friend void h5_write(h5::group g, op_desc const &op) {
      h5_write(g, "block", op.block_index);
      h5_write(g, "inner", op.inner_index);
      h5_write(g, "dagger", op.dagger);
    }

    friend void h5_read(h5::group g, op_desc &op) {
      h5_read(g, "block", op.block_index);
      h5_read(g, "inner", op.inner_index);
      h5_read(g, "dagger", op.dagger);
    }

    bool operator==(op_desc const &op) const { return (block_index == op.block_index && inner_index == op.inner_index
      && dagger == op.dagger && linear_index == op.linear_index);
    }
  };

  // The configuration of the Monte Carlo
  struct configuration {

    bool operator==(configuration const &config) const { return (beta_ == config.beta_ && oplist == config.oplist); }

    // a map associating an operator to an imaginary time
    using oplist_t = std::map<time_pt, op_desc, std::greater<time_pt>>;

#ifdef SAVE_CONFIGS
    configuration(double beta) : beta_(beta), id(0), configs_hfile("configs.h5", exists("configs.h5") ? 'a' : 'w') {
      if (NUM_CONFIGS_TO_SAVE > 0) h5_write(configs_hfile, "c_0", *this);
    }
    ~configuration() { configs_hfile.close(); }
#else
    configuration(double beta) : beta_(beta), id(0) {}
    configuration() : beta_(double(0.)) {}
#endif

    double beta() const { return beta_; }
    int size() const { return oplist.size(); }

    void insert(time_pt tau, op_desc op) { oplist.insert({tau, op}); }
    void replace(time_pt tau, op_desc op) { oplist[tau] = op; }
    void erase(time_pt const &t) { oplist.erase(t); }
    void clear() { oplist.clear(); }

    oplist_t::iterator begin() { return oplist.begin(); }
    oplist_t::iterator end() { return oplist.end(); }
    oplist_t::const_iterator begin() const { return oplist.begin(); }
    oplist_t::const_iterator end() const { return oplist.end(); }

    friend std::ostream &operator<<(std::ostream &out, configuration const &c) {
      for (auto const &op : c) out << "tau = " << op.first << " : " << op.second << std::endl;
      return out;
    }

    static std::string hdf5_format() { return "Configuration"; }

    // Writing of configuration out to a h5 for e.g. plotting
    friend void h5_write(h5::group conf, std::string conf_group_name, configuration const &c) {
      auto beta = c.beta();
      auto t1 = time_pt(1,beta);
      h5::group conf_group = conf.create_group(conf_group_name);
      write_hdf5_format(conf_group, c);
      h5_write(conf_group, "beta", beta);
      for (auto const &op : c) {
        // create group for given tau
        auto tau_group_name        = std::to_string(double(op.first));
        h5::group tau_group        = conf_group.create_group(tau_group_name);
        // in tau subgroup, write operator info
        h5_write(tau_group, op.second);
        h5_write(tau_group, "n", floor_div(op.first, t1));
      }
    }

    friend void h5_read(h5::group conf, std::string conf_group_name, configuration &c) {
      h5::group conf_group = conf.open_group(conf_group_name);
      op_desc op;
      uint64_t n;
      double beta;
      h5_read(conf_group, "beta", beta);
      for (auto &sgrp : conf_group.get_all_subgroup_names()) {
        auto tau_group = conf_group.open_group(sgrp);
        h5_read(tau_group, op);
        h5_read(tau_group, "n", n);
        c.insert(time_pt(n,beta), op);
      }
    }

    long get_id() const { return id; } // Get the id of the current configuration
    void finalize() {
      id++;
#ifdef SAVE_CONFIGS
      if (id < NUM_CONFIGS_TO_SAVE) h5_write(configs_hfile, "c_" + std::to_string(id), *this);
#endif
    }

    private:
    double beta_;
    long id; // configuration id, for debug purposes
    oplist_t oplist;

#ifdef SAVE_CONFIGS
    // HDF5 file to save configurations
    h5::file configs_hfile;
#endif
  };
}
