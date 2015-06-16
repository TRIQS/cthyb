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
#include "impurity_trace.hpp"

namespace cthyb {

struct measure_state_contrib_hist {
 using mc_sign_type = std::complex<double>;

 impurity_trace & imp_tr;
 arrays::vector<double> &contrib;
 double z = 0;

 measure_state_contrib_hist(impurity_trace & imp_tr, arrays::vector<double>& contrib)
    : imp_tr(imp_tr), contrib(contrib) {
  int n_eigstates = imp_tr.n_eigstates;
  contrib.resize(n_eigstates);
  contrib() = 0.0;
 }
 // --------------------

 void accumulate(mc_sign_type s) {
  z += 1;
  imp_tr.estimate(-1,0);
  auto tot = sum(imp_tr.state_contrib);
  if (tot != 0.0)
   this->contrib += imp_tr.state_contrib / tot;
 }
  
 // ---------------------------------------------

 void collect_results(boost::mpi::communicator const& c) {
  arrays::vector<double> contrib_total(imp_tr.n_eigstates);
  std::string s="histo_state_contrib_to_trace.dat";
  for (int i=0; i<contrib.size(); i++) {
   boost::mpi::all_reduce(c, this->contrib[i], contrib_total[i], std::c14::plus<>());
  }
  contrib_total /= z;
  if (c.rank() == 0) {
   std::ofstream f(s);
   size_t i = 0;
   for (auto const& x : contrib_total) {
    f << (i++) << "  " << x << "  " << std::endl;
   }
  }
 }

};
 
}
