#include "solver_core.hpp"
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/gfs.hpp>

using namespace cthyb;
using triqs::utility::many_body_operator;
using triqs::utility::c;
using triqs::utility::c_dag;
using triqs::utility::n;
using namespace triqs::gfs;
using indices_type = triqs::utility::many_body_operator<double>::indices_t;

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
  std::map<std::string, indices_type> gf_struct{{"tot",{0,1}}};

  // Construct CTQMC solver
  solver_core solver(beta, gf_struct, 1025, 2500);

  // Set G0
  triqs::clef::placeholder<0> om_;
  auto g0_iw = gf<imfreq>{{beta, Fermion}, {2,2}};
  g0_iw(om_) << om_ + mu - ( V*V / (om_ - epsilon) + V*V / (om_ + epsilon));
  solver.G0_iw()[0] = triqs::gfs::inverse( g0_iw );

  // Solve parameters
  int n_cycles = 5000;
  auto p = solve_parameters_t(H, n_cycles);
  p.random_name = "";
  p.random_seed = 123 * rank + 567;
  p.max_time = -1;
  p.length_cycle = 50;
  p.n_warmup_cycles = 50;
  p.n_cycles = 5000;
  p.quantum_numbers = qn;
  p.partition_method = "quantum_numbers";
 
  // Solve!
  solver.solve(p);
  
  // Save the results
  if(rank==0){
    triqs::h5::file G_file("anderson_qn.output.h5",H5F_ACC_TRUNC);
    h5_write(G_file,"G_tot",solver.G_tau()[0]);
  }

  return 0;

}
