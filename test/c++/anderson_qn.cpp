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

  // define operators
  auto H = U*n("tot",0)*n("tot",1) + h*n("tot",0) - h*n("tot",1);
  // quantum numbers
  std::vector<many_body_operator<double>> qn;
  qn.push_back(n("tot",0));
  qn.push_back(n("tot",1));
  // gf structure
  std::map<std::string, std::vector<int>> gf_struct{{"tot",{0,1}}};

  // Construct CTQMC solver
  ctqmc solver(beta, gf_struct, 1025, 10001);

  // Set G0
  triqs::clef::placeholder<0> om_;
  auto g0_iw = gf<imfreq>{{beta, Fermion}, {2,2}};
  g0_iw(om_) << om_ + mu - ( V*V / (om_ - epsilon) + V*V / (om_ + epsilon));
  solver.G0_iw_view()[0] = triqs::gfs::inverse( g0_iw );

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
    triqs::h5::file G_file("anderson_qn.output.h5",H5F_ACC_TRUNC);
    h5_write(G_file,"G_tot",solver.G_tau_view()[0]);
  }

  return 0;

}
