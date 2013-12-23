#pragma once
#include <triqs/parameters.hpp>
#include <triqs/mc_tools.hpp>
#include <triqs/utility/callbacks.hpp>
#include <boost/mpi/communicator.hpp>

#include "./qmc_data.hpp"

namespace cthyb_krylov {

using namespace triqs::utility;

class ctqmc_krylov {

 using mc_sign_type = std::complex<double>;

 block_gf<imtime> deltat, gt; // Green's function containers: imaginary-time Green's functions
 boost::mpi::communicator c;  // define the communicator, here MPI_COMM_WORLD
 sorted_spaces sosp;          // Diagonalization of the atomic problem

 public:
 using real_operator_t = many_body_operator<double>;

 ctqmc_krylov(parameters p_in, real_operator_t const& h_loc, std::vector<real_operator_t> const& quantum_numbers,
              fundamental_operator_set const& fops, std::vector<block_desc_t> const& block_structure);

 void solve(parameters p_in);

 // input containers
 block_gf_view<imtime> deltat_view() { return deltat; }

 // imaginary-time measurements
 block_gf_view<imtime> gt_view() { return gt; }

 // specify all required and optional parameters and generate help from them
 parameter_defaults constructor_defaults() const;
 parameter_defaults solve_defaults() const;
 void help() const;
};
}
