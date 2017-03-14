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
#include "qmc_data.hpp"
#include "triqs/statistics/histograms.hpp"

namespace cthyb {

struct measure_perturbation_hist {

 qmc_data const& data;
 int block_index;
 statistics::histogram& histo_perturbation_order;

 measure_perturbation_hist(int block_index, qmc_data const& data, statistics::histogram& hist)
    : data(data), block_index(block_index), histo_perturbation_order(hist) {
  histo_perturbation_order = {0, 1000};
 }

 void accumulate(mc_weight_t s) { histo_perturbation_order << data.dets[block_index].size(); }

 void collect_results(triqs::mpi::communicator const& c) {
  histo_perturbation_order = mpi_all_reduce(histo_perturbation_order, c);
 }
};

// -----------------------------------------------------------------------------
struct measure_perturbation_hist_total {

 qmc_data const& data;
 statistics::histogram& histo_perturbation_order;

 measure_perturbation_hist_total(qmc_data const& data, statistics::histogram& hist) : data(data), histo_perturbation_order(hist) {
  histo_perturbation_order = {0, 1000};
 }

 void accumulate(mc_weight_t s) { histo_perturbation_order << data.config.size() / 2; }

 void collect_results(triqs::mpi::communicator const& c) {
  histo_perturbation_order = mpi_all_reduce(histo_perturbation_order, c);
 }
};
}
