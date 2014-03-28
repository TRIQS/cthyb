#pragma once
#include <triqs/mc_tools.hpp>
#include <triqs/utility/callbacks.hpp>
#include <boost/mpi/communicator.hpp>

#include "./qmc_data.hpp"

namespace cthyb_matrix {

using namespace triqs::utility;
using mc_weight_type = double;

class ctqmc_matrix {

 using mc_sign_type = mc_weight_type;

 block_gf<imtime> deltat, gt; // Green's function containers: imaginary-time Green's functions
 boost::mpi::communicator c;  // define the communicator, here MPI_COMM_WORLD
 sorted_spaces sosp;          // Diagonalization of the atomic problem

 public:
 using real_operator_t = many_body_operator<double>;

 ctqmc_matrix(parameters p_in, real_operator_t const& h_loc, std::vector<real_operator_t> const& quantum_numbers,
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

 block_gf_view<imtime> atomic_gf(double beta) const { return sosp.atomic_gf(beta); }

};
}
