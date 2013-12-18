#pragma once
#include "qmc_data.hpp"
#include <vector>
#include <functional>

namespace cthyb_krylov {

// Measure imaginary time Green's function (one block)
struct measure_path_analysis {

 typedef std::complex<double> mc_sign_type;

 qmc_data const& data;
 triqs::mc_tools::histogram_binned histo;
 std::string stats_file_name;

 // ---------------------------------------------------------

 measure_path_analysis(qmc_data const& data, std::string const& stats_file_name = "path_histo")
    : data(data), stats_file_name(stats_file_name), histo(0, data.config.beta()+ 0.1,100) {}

 // ---------------------------------------------------------

 void accumulate(mc_sign_type s) {

  auto _begin = data.config.oplist.crbegin();
  auto _end = data.config.oplist.crend();
  double beta = data.config.beta();

  // auto last_tau = config.beta();
  int n_blocks = data.sosp.n_subspaces();
  std::vector<double> len(n_blocks, 0);

  std::vector<int> blo(n_blocks);
  for (int u = 0; u < n_blocks; ++u) blo[u] = u;

  // do the first exp
  double dtau = (_begin == _end ? data.config.beta() : double(_begin->first));
  for (int n = 0; n < n_blocks; ++n) len[n] += dtau;

  for (auto it = _begin; it != _end;) { // do nothing if no operator
   auto it1 = it;
   ++it;
   dtau = (it == _end ? data.config.beta() : double(it->first)) - double(it1->first);
   for (int n = 0; n < n_blocks; ++n) {
    if (blo[n] == -1) continue;
    blo[n] = data.sosp.fundamental_operator_connect_from_linear_index(it1->second.dagger, it1->second.linear_index, blo[n]);
    if (blo[n] != -1) len[blo[n]] += dtau;
   }
  }

  //std::cout  << len[n_blocks-1] << std::endl ;
  //if (len[n_blocks-1]>= beta) 
  histo << len[n_blocks -1] ; // CHANGE THIS
   //for (int n = 0; n < n_blocks; ++n) histo << len[n];  
 }

 // ---------------------------------------------------------
 void collect_results(boost::mpi::communicator const& c) { histo.save(std::cout); }
};
}

