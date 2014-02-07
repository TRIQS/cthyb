#pragma once
#include "qmc_data.hpp"
#include <boost/mpi/collectives.hpp>
#include "triqs/statistics/histograms.hpp"

namespace cthyb_krylov {

using namespace triqs::gfs;

// Measure imaginary time Green's function (one block)
struct measure_perturbation_hist {
 using mc_sign_type = std::complex<double>;

 statistics::histogram histo_perturbation_order;
 qmc_data const& data;
 int block_index;

 measure_perturbation_hist(int block_index, qmc_data const& data, std::string hist_file_name)
    : data(data), block_index(block_index), histo_perturbation_order{100, hist_file_name} {
 }
 // --------------------

 void accumulate(mc_sign_type s) {

  histo_perturbation_order << data.dets[block_index].size();
 }
 // ---------------------------------------------

 void collect_results(boost::mpi::communicator const& c) {
 }
};
}
