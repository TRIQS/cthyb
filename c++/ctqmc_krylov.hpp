#ifndef TRIQS_CTQMC_KRYLOV_HPP
#define TRIQS_CTQMC_KRYLOV_HPP

#include <triqs/utility/first_include.hpp>
#include <triqs/parameters.hpp>
#include <triqs/gfs/imtime.hpp>
#include <triqs/gfs/imfreq.hpp>
#include <triqs/gfs/block.hpp>
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <boost/mpi/communicator.hpp>

#include "operator.hpp"
#include "qmc_data.hpp"

#include <vector>

using namespace triqs::gfs;
using namespace triqs::utility;
using triqs::arrays::make_shape;

namespace cthyb_krylov {

class ctqmc_krylov {

 typedef std::complex<double> mc_sign_type;

 // Green's function containers

 // imaginary-time Green's functions
 gf<block_index, gf<imtime>> deltat;
 gf<block_index, gf<imtime>> gt;

 // define the communicator, here MPI_COMM_WORLD
 boost::mpi::communicator c;

 // Information about the atomic problem
 sorted_spaces sosp;

 public:
 using real_operator_t = many_body_operator<double>;

 ctqmc_krylov(parameters p_in, real_operator_t const& h_loc, std::vector<real_operator_t> const& quantum_numbers,
              fundamental_operator_set const& fops, std::vector<block_desc_t> const& block_structure)
    : sosp(h_loc, quantum_numbers, fops, block_structure) {
  p_in.update(constructor_defaults());
  auto const& params = p_in;

  std::vector<std::string> block_names;
  std::vector<gf<imtime>> deltat_blocks;
  std::vector<gf<imtime>> gt_blocks;

  for (auto const& block : block_structure) {
   block_names.push_back(block.name);

   auto shape = make_shape(block.indices.size(), block.indices.size());

   deltat_blocks.push_back(gf<imtime>{{params["beta"], Fermion, params["n_tau_delta"], half_bins}, shape});
   gt_blocks.push_back(gf<imtime>{{params["beta"], Fermion, params["n_tau_g"], half_bins}, shape});
  }

  deltat = make_block_gf(block_names, deltat_blocks);
  gt = make_block_gf(block_names, gt_blocks);
 }

 void solve(parameters p_in);

 // input containers
 gf_view<block_index, gf<imtime>> deltat_view() { return deltat; }

 // imaginary-time measurements
 gf_view<block_index, gf<imtime>> gt_view() { return gt; }

 // specify all required and optional parameters and generate help from them
 parameter_defaults constructor_defaults() const;
 parameter_defaults solve_defaults() const;
 void help() const;
};
}
#endif
