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
#include <triqs/utility/tau_t.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/atom_diag/functions.hpp>
#include <triqs/utility/macros.hpp>

#include <h5/h5.hpp>

#include <map>

namespace triqs_cthyb {

  using triqs::utility::tau_t;

  /// Description of a creation/annihilation operator.
  struct op_desc {
    /// Block index of the operator.
    int block_index;

    /// Inner index within the block.
    int inner_index;

    /// Whether the operator is a dagger (creation operator).
    bool dagger;

    /// Cumulative (linear) index.
    long linear_index;

    friend std::ostream &operator<<(std::ostream &out, op_desc const &op) {
      out << (op.dagger ? "Cdag(" : "C(") << op.block_index << "," << op.inner_index << ")";
      return out;
    }

    static std::string hdf5_format() { return "op_desc"; }

    friend void h5_write(h5::group g, std::string const &name, op_desc const &op) {
      auto gr = g.create_group(name);
      h5::write_hdf5_format(gr, op); // NOLINT (slicing is intended)
      h5::write(gr, "block", op.block_index);
      h5::write(gr, "inner", op.inner_index);
      h5::write(gr, "dagger", op.dagger);
      h5::write(gr, "linear_index", op.linear_index);
    }

    friend void h5_read(h5::group g, std::string const &name, op_desc &op) {
      h5::group gr = g.open_group(name);
      h5::assert_hdf5_format(gr, op);
      h5::read(g, "block", op.block_index);
      h5::read(g, "inner", op.inner_index);
      h5::read(g, "dagger", op.dagger);
      h5::read(g, "linear_index", op.linear_index);
    }

    bool operator==(op_desc const &op) const = default;
  };

  /// Configuration of the Monte Carlo simulation (operators on the imaginary-time line).
  struct configuration {

    bool operator==(configuration const &config) const { return (beta_ == config.beta_ && oplist_ == config.oplist_); }

    // a map associating an operator to an imaginary time
    using oplist_t = std::map<tau_t, op_desc, std::greater<tau_t>>;

#ifdef SAVE_CONFIGS
    configuration(double beta, long id = 0, oplist_t oplist = {})
       : beta_(beta), id_(id), oplist_(oplist), configs_hfile("configs.h5", exists("configs.h5") ? 'a' : 'w') {
      if (NUM_CONFIGS_TO_SAVE > 0) h5_write(configs_hfile, "c_0", *this);
    }
    ~configuration() { configs_hfile.close(); }
#else
    configuration(double beta, long id = 0) : beta_(beta), id_(id) {}
    C2PY_IGNORE configuration(double beta, long id, oplist_t oplist) : beta_(beta), id_(id), oplist_(oplist) {}
#endif

    /// Inverse temperature \f$ \beta \f$.
    C2PY_PROPERTY_GET(beta) double beta() const { return beta_; }
    auto size() const { return oplist_.size(); }

    /**
     * @brief Insert a given operator at a given imaginary time.
     * 
     * @param tau Imaginary time at which to insert the operator.
     * @param op Description of the operator to insert.
     */
    void insert(tau_t tau, op_desc op) { oplist_.insert({tau, op}); }

    /**
     * @brief Replace an existing operator at a given imaginary time with a new one.
     * 
     * @param tau Imaginary time at which to replace the operator.
     * @param op Description of the operator to insert.
     */
    void replace(tau_t tau, op_desc op) { oplist_[tau] = op; }

    /**
     * @brief Erase the operator at a given imaginary time.
     * @param tau Imaginary time at which to erase the operator.
     */
    void erase(tau_t const &t) { oplist_.erase(t); }

    /// Clear the configuration (remove all operators).
    void clear() { oplist_.clear(); }

    oplist_t::iterator begin() { return oplist_.begin(); }
    oplist_t::iterator end() { return oplist_.end(); }
    oplist_t::const_iterator begin() const { return oplist_.begin(); }
    oplist_t::const_iterator end() const { return oplist_.end(); }

    friend std::ostream &operator<<(std::ostream &out, configuration const &c) {
      for (auto const &op : c) out << "tau = " << op.first << " : " << op.second << std::endl;
      return out;
    }

    /// HDF5 format string for configuration.
    static std::string hdf5_format() { return "CTHYB_Configuration"; }

    /// Write a configuration to an hdf5 file.
    friend void h5_write(h5::group g, std::string const &name, configuration const &c) {
      h5::group gr = g.create_group(name);
      h5::write_hdf5_format(gr, c); // NOLINT (slicing is intended)
      h5::write(gr, "beta", c.beta_);
      h5::write(gr, "id", c.id_);
      h5::write(gr, "oplist", c.oplist_);
    }

    /// Read a configuration from an hdf5 file.
    C2PY_IGNORE static configuration h5_read_construct(h5::group g, std::string const &name) {
      h5::group gr = g.open_group(name);
      h5::assert_hdf5_format<configuration>(gr);
      auto beta   = h5::read<double>(gr, "beta");
      auto id     = h5::read<long>(gr, "id");
      auto oplist = h5::read<oplist_t>(gr, "oplist");
      return configuration(beta, id, std::move(oplist));
    }

    /// Get the ID of the current configuration (for debug purposes).
    long get_id() const { return id_; } // Get the id of the current configuration

    /// Finalize the configuration after a Monte Carlo move (increment the ID and save the configuration if needed).
    void finalize() {
      id_++;
#ifdef SAVE_CONFIGS
      if (id < NUM_CONFIGS_TO_SAVE) h5_write(configs_hfile, "c_" + std::to_string(id), *this);
#endif
    }

    private:
    double beta_;
    long id_; // configuration id, for debug purposes
    oplist_t oplist_;

#ifdef SAVE_CONFIGS
    // HDF5 file to save configurations
    h5::file configs_hfile;
#endif
  };
} // namespace triqs_cthyb
