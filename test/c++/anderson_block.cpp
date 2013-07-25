#include "./ctqmc.hpp"
#include "./operator.hpp"
#include "./fundamental_operator_set.hpp"
#include <triqs/gf/local/fourier_matsubara.hpp>
#include <triqs/parameters.hpp>
#include <triqs/gf/block.hpp>
#include <triqs/gf/imtime.hpp>
#include <triqs/gf/imfreq.hpp>

using namespace triqs::app::impurity_solvers::ctqmc_krylov;
using triqs::utility::many_body_operator;
using triqs::utility::c;
using triqs::utility::c_dag;
using triqs::utility::n;
using triqs::utility::parameters;
using triqs::gf::gf;
using triqs::gf::block_index;
using triqs::gf::imfreq;
using triqs::gf::imtime;
using triqs::gf::make_gf;
using triqs::gf::Fermion;

#include <boost/lexical_cast.hpp>

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

  // Put in the class
  parameters p;
  p["beta"] = beta;
  p["max_time"] = -1;
  p["Random_Generator_Name"] = "";
  p["Random_Seed"] = 123 * rank + 567;
  p["Verbosity"] = 3;
  p["Length_Cycle"] = 50;
  p["N_Warmup_Cycles"] = 10;
  p["N_Cycles"] = 5000;
  
  // define operators
  auto H = U*n("up")*n("down") + (-mu+h)*n("up") + (-mu-h)*n("down");

  // quantum numbers
  std::vector<many_body_operator<double,const char*>> qn;

  // basis of operators to use
  fundamental_operator_set<const char *> fops;
  fops.add_operator("up");
  fops.add_operator("down");

  // map indices --> pair of ints
  std::map<std::tuple<const char *>, std::pair<int,int>> my_map;

  // decide wether g has a block structure
  int dim_block, n_blocks;
  std::vector<std::string> block_names;
  my_map[std::make_tuple("up")] = std::make_pair(0,0);
  my_map[std::make_tuple("down")] = std::make_pair(1,0);
  block_names.push_back("up");
  block_names.push_back("down");
  dim_block = 1;
  n_blocks = 2;

  // Green's functions
  auto sha = triqs::arrays::make_shape(dim_block,dim_block);
  auto Delta = make_gf<block_index, gf<imtime>>(block_names, make_gf<imtime>(beta, Fermion, sha, 1000) );
  auto G = make_gf<block_index, gf<imtime>>(block_names, make_gf<imtime>(beta, Fermion, sha, 1000) );

  // Set hybridization function
  triqs::clef::placeholder<0> om_;
  auto delta_w = make_gf<imfreq>(beta, Fermion, sha);
  delta_w(om_) << V*V / (om_ - epsilon) + V*V / (om_ + epsilon);
  for (int i=0; i<n_blocks; i++) Delta()[i] = triqs::gf::lazy_inverse_fourier(delta_w);

  // Construct CTQMC solver
  ctqmc krylov_ctqmc(p, H, qn, fops, my_map, G, Delta);
  
  // Solve!
  krylov_ctqmc.solve();
  
  // Save the results
  if(rank==0){
    H5::H5File G_file("anderson_block.output.h5",H5F_ACC_TRUNC);
    h5_write(G_file,"G_up",G[0]);
    h5_write(G_file,"G_down",G[1]);
  }

  return 0;

}
