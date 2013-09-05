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
using triqs::gfs::imfreq;
using triqs::gfs::imfreq;
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
  p["random_name"] = "";
  p["random_seed"] = 123 * rank + 567;
  p["verbosity"] = 3;
  p["length_cycle"] = 50;
  p["n_warmup_cycles"] = 50;
  p["n_cycles"] = 500;
  p["n_tau_delta"] = 1000;
  p["n_tau_g"] = 1000;
  p["krylov_bs_use_cutoff"] = true;
  p["krylov_bs_prob_cutoff"] = .0;

  // basis of operators to use  
  fundamental_operator_set<const char *,int> fops;
  for(int o = 0; o < num_orbitals; ++o){
      fops.add_operator("up",o);
      fops.add_operator("down",o);
  }
  
  // block structure of GF
  std::vector<block_desc_t<const char *,int>> block_structure;
  block_structure.push_back({"up",{}});
  block_structure.push_back({"down",{}});
  for(int o = 0; o < num_orbitals; ++o){
      block_structure[0].indices.push_back(std::make_tuple("up",o));
      block_structure[1].indices.push_back(std::make_tuple("down",o));
  }

  // Hamiltonian
  many_body_operator<double,const char*,int> H;
  for(int o = 0; o < num_orbitals; ++o){
      H += -mu*(n("up",o) + n("down",o));
  }
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

  // quantum numbers
  std::vector<many_body_operator<double,const char*,int>> qn;
  qn.resize(2);
  for(int o = 0; o < num_orbitals; ++o){
    qn[0] += n("up",o);
    qn[1] += n("down",o);
  }

  // Construct CTQMC solver
  ctqmc_krylov solver(p, H, qn, fops, block_structure);
  
  // Set hybridization function
  auto delta_w = make_gf<imfreq>(beta, Fermion, make_shape(num_orbitals,num_orbitals));
  
  triqs::clef::placeholder<0> om_;
  auto term = make_gf<imfreq>(beta, Fermion, make_shape(num_orbitals,num_orbitals));  
  for(int j=0; j < epsilon.size(); ++j){
      term(om_) << 1.0/(om_ - epsilon[j]);
      
      matrix<std::complex<double>> m(num_orbitals,num_orbitals);
      for(std::size_t w_index = 0; w_index < term.mesh().size(); ++w_index){
          m = term.data()(w_index,ellipsis());
          m = _conj(V[j]) * m * V[j];
          term.data()(w_index,ellipsis()) = m;
      }
      for(int tail_o = term.singularity().order_min();
              tail_o <= term.singularity().order_max(); ++tail_o){
          m = term.singularity()(tail_o);
          term.singularity()(tail_o) = _conj(V[j]) * m * V[j];
      }
      delta_w = delta_w + term;
  }
  
  solver.deltat_view()[0] = triqs::gfs::lazy_inverse_fourier(delta_w);
  solver.deltat_view()[1] = triqs::gfs::lazy_inverse_fourier(delta_w);
  
  // Solve!
  solver.solve(p);
  
  // Save the results
  if(rank==0){
    H5::H5File G_file("kanamori_offdiag_qn.output.h5",H5F_ACC_TRUNC);
    h5_write(G_file,"G_up",solver.gt_view()[0]);
    h5_write(G_file,"G_down",solver.gt_view()[1]);
  }

  return 0;

}
