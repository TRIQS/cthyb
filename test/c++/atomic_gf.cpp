#include "solver_core.hpp"
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/draft/hilbert_space_tools/fundamental_operator_set.hpp>
#include <triqs/parameters.hpp>
#include <triqs/gfs/imfreq.hpp>
#include <triqs/gfs/block.hpp>

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
  double U0 = 2.0;
  double U1 = 2.0;
  double ed0 = -1.1;
  double ed1 = -0.9;
  double V = 1.5;

  // define operators
  auto H = U0*n("up",0)*n("dn",0) + U1*n("up",1)*n("dn",1);
  H += ed0*(n("up",0) + n("dn",0)) + ed1*(n("up",1) + n("dn",1));
  H += V*(c_dag("up",0)*c("up",1) + c_dag("up",1)*c("up",0) + c_dag("dn",0)*c("dn",1) + c_dag("dn",1)*c("dn",0));

  std::map<std::string, std::vector<int>> gf_struct{{"up",{0,1}},{"dn",{0,1}}};

  // Construct CTQMC solver
  solver_core solver(beta, gf_struct, 1025, 2051);

  // Solve parameters
  auto p = solver_core::solve_parameters();
  p["length_cycle"] = 1;
  p["n_warmup_cycles"] = 1;
  p["n_cycles"] = 1;

  triqs::clef::placeholder<0> om_;
  auto g0_iw = gf<imfreq>{{beta, Fermion}, {2,2}};
  g0_iw(om_) << om_ + 0.0;
  solver.G0_iw()[0] = triqs::gfs::inverse(g0_iw);
  solver.G0_iw()[1] = triqs::gfs::inverse(g0_iw);

  // Solve!
  solver.solve(H, p);

  // Save the results
  if(rank==0){
    triqs::h5::file G_file("atomic_gf.output.h5",H5F_ACC_TRUNC);
    h5_write(G_file,"G_tau",solver.atomic_gf());
  }

  return 0;
}
