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
using triqs::gf::imfreq;
using triqs::gf::imfreq;
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
  int num_orbitals = 2;
  double mu = 1.0;
  double U = 2.0;
  double J = 0.2;
  bool n_n_only = false;
  
  // Poles of delta
  std::vector<double> epsilon;
  epsilon.push_back(2.3);
  epsilon.push_back(-2.3);
  
  // Hybridization matrices
  std::vector<triqs::arrays::matrix<std::complex<double>>> V;
  V.emplace_back(num_orbitals,num_orbitals);
  V.emplace_back(num_orbitals,num_orbitals);
  
  for(int o1 = 0; o1 < num_orbitals; ++o1)
  for(int o2 = 0; o2 < num_orbitals; ++o2){
      V[0](o1,o2) = (o1 == o2 ? 1.0 : 0.1);
      V[1](o1,o2) = (o1 == o2 ? 1.0 : 0.1);
  }
    
  assert(epsilon.size() == V.size());
  
  // Put in the class
  parameters p;
  p["beta"] = beta;
  p["max_time"] = -1;
  p["Random_Generator_Name"] = "";
  p["Random_Seed"] = 123 * rank + 567;
  p["Verbosity"] = 3;
  p["Length_Cycle"] = 50;
  p["N_Warmup_Cycles"] = 50;
  p["N_Cycles"] = 1000;

  // basis of operators to use  
  fundamental_operator_set<int, const char *> fops;
  for(int o = 0; o < num_orbitals; ++o){
      fops.add_operator(o,"up");
      fops.add_operator(o,"down");
  }

  // Hamiltonian
  many_body_operator<double,int, const char*> H;
  for(int o = 0; o < num_orbitals; ++o){
      H += -mu*(n(o,"up") + n(o,"down"));
  }
  for(int o = 0; o < num_orbitals; ++o){
      H += U *n(o,"up")*n(o,"down");
  }
  for(int o1 = 0; o1 < num_orbitals; ++o1)
  for(int o2 = 0; o2 < num_orbitals; ++o2){
      if(o1==o2) continue;
      H += (U-2*J)*n(o1,"up")*n(o2,"down");
  }
  for(int o1 = 0; o1 < num_orbitals; ++o1)
  for(int o2 = 0; o2 < num_orbitals; ++o2){
      if(o2>=o1) continue;
      H += (U-3*J)*n(o1,"up")*n(o2,"up");
      H += (U-3*J)*n(o1,"down")*n(o2,"down");
  }
  
  if(!n_n_only){ // spin flips and pair hopping
    for(int o1 = 0; o1 < num_orbitals; ++o1)
    for(int o2 = 0; o2 < num_orbitals; ++o2){
        if(o1==o2) continue;
        H += -J*c_dag(o1,"up")*c_dag(o1,"down")*c(o2,"up")*c(o2,"down");
        H += -J*c_dag(o1,"up")*c_dag(o2,"down")*c(o2,"up")*c(o1,"down");
    }
  }

  // quantum numbers
  std::vector<many_body_operator<double,int, const char*>> qn;
  qn.resize(2);
  for(int o = 0; o < num_orbitals; ++o){
    qn[0] += n(o,"up");
    qn[1] += n(o,"down");
  }
  
  // map indices --> pair of ints
  std::map<std::tuple<int,const char *>, std::pair<int,int>> my_map;
  for(int o = 0; o < num_orbitals; ++o){
      my_map[std::make_tuple(o,"up")] = std::make_pair(0,o);
      my_map[std::make_tuple(o,"down")] = std::make_pair(1,o);
  }

  // Green's functions
  std::vector<std::string> block_names;
  block_names.push_back("up");
  block_names.push_back("down");
  auto sha1 = triqs::arrays::make_shape(num_orbitals,num_orbitals);

  auto Delta = make_gf<block_index, gf<imtime>>(block_names, make_gf<imtime>(beta, Fermion, sha1, 1000) );
  auto G = make_gf<block_index, gf<imtime>>(block_names, make_gf<imtime>(beta, Fermion, sha1, 1000) );

  // Set hybridization function
  auto delta_w = make_gf<imfreq>(beta, Fermion, sha1);
    
  auto w_mesh = delta_w.mesh();
  for(std::size_t w_index = 0; w_index < w_mesh.size(); ++w_index){
      auto iw = w_mesh.index_to_point(w_index);
      
      auto m = delta_w(w_index);
      for(int j=0; j < epsilon.size(); ++j){
          for(int o = 0; o < num_orbitals; ++o) m(o,o) = 1.0/(iw - epsilon[j]);
          m = _conj(V[j]) * m * V[j];
      }
  }
  
  Delta()[0] = triqs::gf::lazy_inverse_fourier(delta_w);
  Delta()[1] = triqs::gf::lazy_inverse_fourier(delta_w);

  // Construct CTQMC solver
  ctqmc krylov_ctqmc(p, H, qn, fops, my_map, G, Delta);
  
  // Solve!
  krylov_ctqmc.solve();
  
  // Save the results
  if(rank==0){
    H5::H5File G_file("kanamori_offdiag_qn.output.h5",H5F_ACC_TRUNC);
    h5_write(G_file,"G_up",G[0]);
    h5_write(G_file,"G_down",G[1]);
  }

  return 0;

}
