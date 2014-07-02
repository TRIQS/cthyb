#include "ctqmc.hpp"
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/draft/hilbert_space_tools/fundamental_operator_set.hpp>
#include <triqs/gfs/local/fourier_matsubara.hpp>
#include <triqs/parameters.hpp>
#include <triqs/gfs/block.hpp>
#include <triqs/gfs/imtime.hpp>
#include <triqs/gfs/imfreq.hpp>

using namespace cthyb;
using triqs::utility::many_body_operator;
using triqs::utility::c;
using triqs::utility::c_dag;
using triqs::utility::n;
using triqs::params::parameters;
using namespace triqs::gfs;

int main(int argc, char* argv[]) {

  std::cout << "Welcome to the CTHYB solver\n";

  // Initialize mpi
  boost::mpi::environment env(argc, argv);
  int rank;
  {
    boost::mpi::communicator c;
    rank = c.rank();
  }

  // Parameters
  double beta = 10.0;
  double U = 2.0;
  double mu = 1.0;
  double h = 0.0;
  double V = 1.0;
  double epsilon = 2.3;

  // define operators and QN
  auto H = U*n("up",0)*n("down",0) + (-mu+h)*n("up",0) + (-mu-h)*n("down",0);
  std::vector<many_body_operator<double>> qn;
  std::map<std::string, std::vector<int>> gf_struct{{"up",{0}},{"down",{0}}};

  // Construct CTQMC solver
  ctqmc solver(beta, gf_struct, 1000, 1000);

  // Set hybridization function
  triqs::clef::placeholder<0> om_;
  auto delta_iw = gf<imfreq>{{beta, Fermion}, {1,1}};
  delta_iw(om_) << V*V / (om_ - epsilon) + V*V / (om_ + epsilon);  
  for (int bl=0; bl<2; ++bl) solver.Delta_tau_view()[bl] = triqs::gfs::inverse_fourier(delta_iw);

  // Solve parameters
  auto p = ctqmc::solve_parameters();
  p["random_name"] = "";
  p["random_seed"] = 123 * rank + 567;
  p["max_time"] = -1;
  p["verbosity"] = 3;
  p["length_cycle"] = 50;
  p["n_warmup_cycles"] = 10;
  p["n_cycles"] = 5000;

  // Solve!
  solver.solve(H, p, qn, true);
  
  // Save the results
  if(rank==0){
    triqs::h5::file G_file("anderson_block.output.h5",H5F_ACC_TRUNC);
    h5_write(G_file,"G_up",solver.G_tau_view()[0]);
    h5_write(G_file,"G_down",solver.G_tau_view()[1]);
  }

  return 0;

}
