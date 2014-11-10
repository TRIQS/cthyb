#include "solver_core.hpp"
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/draft/hilbert_space_tools/fundamental_operator_set.hpp>
#include <triqs/gfs.hpp>

using namespace cthyb;
using triqs::utility::many_body_operator;
using triqs::utility::c;
using triqs::utility::c_dag;
using triqs::utility::n;
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
  
  // Hamiltonian
  many_body_operator<double> H;
  for(int o = 0; o < num_orbitals; ++o){
      H += U *n("up",o)*n("down",o);
  }
  for(int o1 = 0; o1 < num_orbitals; ++o1)
  for(int o2 = 0; o2 < num_orbitals; ++o2){
      if(o1==o2) continue;
      H += (U-2*J)*n("up",o1)*n("down",o2);
  }
  for(int o1 = 0; o1 < num_orbitals; ++o1)
  for(int o2 = 0; o2 < num_orbitals; ++o2){
      if(o2>=o1) continue;
      H += (U-3*J)*n("up",o1)*n("up",o2);
      H += (U-3*J)*n("down",o1)*n("down",o2);
  }
  
  if(!n_n_only){ // spin flips and pair hopping
    for(int o1 = 0; o1 < num_orbitals; ++o1)
    for(int o2 = 0; o2 < num_orbitals; ++o2){
        if(o1==o2) continue;
        H += -J*c_dag("up",o1)*c_dag("down",o1)*c("up",o2)*c("down",o2);
        H += -J*c_dag("up",o1)*c_dag("down",o2)*c("up",o2)*c("down",o1);
    }
  }

  // gf structure
  std::map<std::string, std::vector<int>> gf_struct{{"up",std::vector<int>{}},{"down",std::vector<int>{}}}; 
  for(int o = 0; o < num_orbitals; ++o){
    gf_struct["up"].push_back(o); 
    gf_struct["down"].push_back(o); 
  }

  // Construct CTQMC solver
  solver_core solver(beta, gf_struct, 1025, 2500);

  // Set G0
  auto delta_iw = gf<imfreq>{{beta, Fermion}, {num_orbitals,num_orbitals}};

  auto w_mesh = delta_iw.mesh();
  for(std::size_t w_index = 0; w_index < w_mesh.size(); ++w_index){
      auto iw = w_mesh.index_to_point(w_index);
      
      auto m = delta_iw[w_index];
      for(int j=0; j < epsilon.size(); ++j){
          for(int o = 0; o < num_orbitals; ++o) m(o,o) = 1.0/(iw - epsilon[j]);
          m = conj(V[j]) * m * V[j];
      }
  }

  triqs::clef::placeholder<0> om_;
  auto g0_iw = gf<imfreq>{{beta, Fermion}, {num_orbitals,num_orbitals}};
  g0_iw(om_) << om_ + mu - delta_iw(om_);
  solver.G0_iw()[0] = triqs::gfs::inverse( g0_iw );
  solver.G0_iw()[1] = triqs::gfs::inverse( g0_iw );

  // Solve parameters
  auto n_cycles = 5000;
  auto p = solve_parameters_t(H, n_cycles);
  p.max_time = -1;
  p.random_name = "";
  p.random_seed = 123 * rank + 567;
  p.verbosity = 3;
  p.length_cycle = 50;
  p.n_warmup_cycles = 50;

  // Solve!
  solver.solve(p);
  
  // Save the results
  if(rank==0){
    triqs::h5::file G_file("kanamori_offdiag.output.h5",H5F_ACC_TRUNC);
    h5_write(G_file,"G_up",solver.G_tau()[0]);
    h5_write(G_file,"G_down",solver.G_tau()[1]);
  }

  return 0;

}
