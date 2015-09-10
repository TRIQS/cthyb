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
#include "./measure_density_matrix.hpp"

namespace cthyb {

measure_density_matrix::measure_density_matrix(qmc_data const& data, std::vector<matrix<double>>& density_matrix)
   : data(data), block_dm(density_matrix) {
 block_dm.resize(data.imp_trace.get_density_matrix().size());
 for (int i = 0; i < block_dm.size(); ++i) {
  block_dm[i] = data.imp_trace.get_density_matrix()[i].mat;
  block_dm[i]() = 0;
 }
}
// --------------------

void measure_density_matrix::accumulate(mc_sign_type s) {
 // we assume here that we are in "Norm" mode, i.e. qmc weight is norm, not trace

 // We need to recompute since the density_matrix in the trace is changed at each computatation,
 // in particular at the last failed attempt.
 // So we need to compute it, without any Yee threshold.
 data.imp_trace.compute();
 z += s * data.atomic_reweighting;
 s /= data.atomic_weight; // accumulate matrix / norm since weight is norm * det

 // Careful: there is no reweighting factor here!
 int size = block_dm.size();
 for (int i = 0; i < size; ++i)
  if (data.imp_trace.get_density_matrix()[i].is_valid) {
   block_dm[i] += s * data.imp_trace.get_density_matrix()[i].mat;
  }

// FIXME

 /* // debug trace
 double tr = 0;
 for (int i = 0; i < size; ++i)
  if (data.imp_trace.get_density_matrix()[i].is_valid) {
   tr +=  trace(data.imp_trace.get_density_matrix()[i].mat);
  }

 if (std::abs(tr-data.atomic_reweighting *  data.atomic_weight)> 1.e-15) TRIQS_RUNTIME_ERROR << "Trace of the density matrix acc
 error "<< tr <<"  "<< data.atomic_reweighting *  data.atomic_weight;
 */
}

// ---------------------------------------------

void measure_density_matrix::collect_results(boost::mpi::communicator const& c) {
 auto block_dm2 = block_dm;
 mc_sign_type total_z;

 // boost::mpi::all_reduce(c, block_dm, block_dm2, std::c14::plus<>());
 // mc_sign_type total_z = mpi_all_reduce(z, c);
 boost::mpi::all_reduce(c, z, total_z, std::c14::plus<>());
 for (int i = 0; i < block_dm2.size(); ++i) boost::mpi::all_reduce(c, block_dm2[i], block_dm[i], std::c14::plus<>());

 for (auto& b : block_dm) b = b / real(total_z);

 if (c.rank() != 0) return;

 // Check: the trace of the density matrix must be 1 by construction
 double tr = 0;
 for (auto& B : block_dm) tr += trace(B);
 if (std::abs(tr - 1) > 0.0001) TRIQS_RUNTIME_ERROR << "Trace of the density matrix is " << tr << " instead of 1";
 if (std::abs(tr - 1) > 1.e-13) std::cerr << "Warning :: Trace of the density matrix is " << tr << " instead of 1";

// FIXME
// h5_write( h5::file("density_matrix.h5", 'w'), "density_matrix", block_dm);
}
}
