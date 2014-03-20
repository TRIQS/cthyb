/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011-13 by O. Parcollet
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
#include <map>
#include <fstream>
#include <boost/mpi.hpp>

namespace triqs {
namespace statistics {

 /**
   Histograms for a MonteCarlo run.
   For integer in [0,Size]
   */
 class histogram {
  uint64_t del;
  std::vector<uint64_t> data;
  std::string _dump_file_name;
  int imin;

  public:
  /**
   * Constructs an histogram for int from i_min, to i_max
  */
  histogram(int i_min, int i_max, std::string dump_file_name = "") : del(0), data(i_max - i_min + 1, 0), imin(i_min) {
   activate_dumpfile(dump_file_name);
  }

  histogram(int i_max, std::string dump_file_name = "") : histogram(0, i_max, dump_file_name) {}

  histogram() : histogram(0, 10, "") {}

  ///
  ~histogram() {
   if (!_dump_file_name.empty()) dump(_dump_file_name);
  }

  /// Activate a dump file if string is not "" : the histogram will be saved automatically there at destruction
  void activate_dumpfile(std::string dump_file_name) { _dump_file_name = dump_file_name; }

  /// Accumulate an integer into the histogram
  histogram& operator<<(int i) {
   auto n = i - imin;
   if ((n < 0) || (n >= data.size()))
    ++del;
   else
    ++data[n];
   return *this;
  }

  /// Number of bins
  size_t n_bins() const { return data.size(); }

  uint64_t n_lost_pts() const { return del; }

  ///
  void clear() {
   del = 0;
   for (auto& x : data) x = 0;
  }

  /// Return the normalized histogram (sum = 1)
  std::vector<double> normalize() const {
   std::vector<double> r(data.size());
   double norm = 0;
   for (auto x : data) norm += x;
   size_t i = 0;
   for (auto const& x : data) r[i++] = x / norm;
   return r;
  }

  /** */
  void dump(std::string s) {
   boost::mpi::communicator world;
   if (world.rank() == 0) {
    std::ofstream f(s);
    size_t i = 0;
    double cum=0;
    for (auto const& x : normalize()) {
     cum += x;
     //f << std::setw(4) << i++ << "  " << std::setw(10) << x << " " << std::setw(10) << cum << std::endl;
     f <<  i++ << "  " <<  x << " " <<  cum << std::endl;
    }
    if (del) std::cerr << "Histogram : " << del << " points have been lost !" << std::endl;
   } else {
    //std::cout << "not dumping histo " << s <<" on node " << world.rank() << std::endl;
   }
  }
 };

 //-------------------------------------------------------------------------------
 /**
   Histogram binning on a segment
   */
 class histogram_segment_bin {
  histogram histo;
  double _a, _b, n_bin_over_len;
  std::string _dump_file_name;

  public:
  /** */
  histogram_segment_bin() : histogram_segment_bin(0, 1) {}

  histogram_segment_bin(double a, double b, int n_bins = 1000, std::string dump_file_name = "")
     : histo(n_bins+1), _a(a), _b(b), n_bin_over_len(n_bins / (b - a)) {
   activate_dumpfile(dump_file_name);
   if (_a >= _b) TRIQS_RUNTIME_ERROR << "histogram_segment_bin construction : one must have a<b";
  }

  // TO FIX DUMPING EMPTY.... DOS NOT WORK !
  histogram_segment_bin(histogram_segment_bin const&) = default;
  histogram_segment_bin(histogram_segment_bin&&) = default;
  histogram_segment_bin& operator=(histogram_segment_bin const&) = default;
  histogram_segment_bin& operator=(histogram_segment_bin&&) = default;

  ///
  ~histogram_segment_bin() {
   if (!_dump_file_name.empty()) dump(_dump_file_name);
  }

  /// Bins a double into the histogram
  histogram_segment_bin& operator<<(double x) {
   histo << int(floor((x - _a) * n_bin_over_len));
   return *this;
  }

  /// Activate a dump file if string is not "" : the histogram will be saved automatically there at destruction
  void activate_dumpfile(std::string dump_file_name) { _dump_file_name = dump_file_name; }

  /// Dump into text file
  void dump(std::string s) {
   boost::mpi::communicator world;
   if (world.rank() == 0) {
    std::cout << "DUMPING HISTO" << std::endl;
    std::ofstream f(s);
    size_t i = 0;
    double cum=0;
    for (auto const& x : histo.normalize()) {
     cum +=x;
     f << _a + (i++) / n_bin_over_len << "  " << x << "  " << cum << std::endl;
    }
    if (histo.n_lost_pts() != 0) std::cerr << "Histogram : " << histo.n_lost_pts() << " points have been lost !" << std::endl;
   } else {
    //std::cout << "not dumping histo "<< s << " on node " << world.rank() << std::endl;
   }
  }
 };
}
}

