#include "./ctqmc_krylov.hpp"
#include "./operator.hpp"
#include "./fundamental_operator_set.hpp"
#include <triqs/gfs/local/fourier_matsubara.hpp>
#include <triqs/parameters.hpp>
#include <triqs/gfs/block.hpp>
#include <triqs/gfs/imtime.hpp>
#include <triqs/gfs/imfreq.hpp>

using namespace triqs::app::impurity_solvers::ctqmc_krylov;
using triqs::utility::many_body_operator;
using triqs::utility::c;
using triqs::utility::c_dag;
using triqs::utility::n;
using triqs::utility::parameters;
using triqs::gfs::gf;
using triqs::gfs::block_index;
using triqs::gfs::imfreq;
using triqs::gfs::imtime;
using triqs::gfs::make_gf;
using triqs::gfs::Fermion;

int main(int argc, char* argv[]) {

  std::cout << "Welcome to the Krylov solver\n";

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

  parameters p;
  p["beta"] = beta;
  p["random_name"] = "";
  p["random_seed"] = 123 * rank + 567;
  p["max_time"] = -1;
  p["verbosity"] = 3;
  p["Length_Cycle"] = 50;
  p["n_warmup_cycles"] = 10;
  p["n_cycles"] = 5000;
  p["n_tau_delta"] = 1000;
  p["n_tau_g"] = 1000;
  p["krylov_bs_use_cutoff"] = true;
  p["krylov_bs_prob_cutoff"] = .0;
  
  // define operators
  auto H = U*n("up")*n("down") + (-mu+h)*n("up") + (-mu-h)*n("down");

  // quantum numbers
  std::vector<many_body_operator<double,std::string>> qn;

  // basis of operators to use
  fundamental_operator_set<std::string> fops;
  fops.add_operator("up");
  fops.add_operator("down");
 
  // block structure of GF
  std::vector<block_desc_t<std::string>> block_structure;
  block_structure.push_back({"up",{std::make_tuple("up")}});
  block_structure.push_back({"down",{std::make_tuple("down")}});
  
  // Construct CTQMC solver
  ctqmc_krylov solver(p, H, qn, fops, block_structure);

  // Set hybridization function
  triqs::clef::placeholder<0> om_;
  auto delta_w = make_gf<imfreq>(beta, Fermion, make_shape(1,1));
  delta_w(om_) << V*V / (om_ - epsilon) + V*V / (om_ + epsilon);  
  for (int bl=0; bl<2; ++bl) solver.deltat_view()[bl] = triqs::gfs::lazy_inverse_fourier(delta_w);
  
  // Solve!
  solver.solve(p);
  
  // Save the results
  if(rank==0){
    H5::H5File G_file("anderson_block.output.h5",H5F_ACC_TRUNC);
    h5_write(G_file,"G_up",solver.gt_view()[0]);
    h5_write(G_file,"G_down",solver.gt_view()[1]);
  }
  
  return 0;

}
