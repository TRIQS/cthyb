#include <triqs/utility/first_include.hpp>
#include <vector>
#include <iostream>
#include <complex>

#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>
#include <triqs/arrays/asserts.hpp>

using triqs::arrays::vector;
using triqs::arrays::matrix;
using triqs::arrays::dotc;

std::size_t get_space_dim(vector<double> const& st)
{
    return st.size();
}

vector<double> make_zero_state(vector<double> const& st)
{
    vector<double> zero_st(st.size());
    zero_st() = 0;
    return zero_st;
}

#include "./krylov_worker.hpp"
#include "./exp_h_worker.hpp"

using namespace triqs::app::impurity_solvers::ctqmc_krylov;


int main() {

    // test krylov_worker
    matrix<double> h(5,5);  // Hamiltonian matrix
    h(0,0) = 1;
    h(1,1) = 2;
    h(2,2) = 3;
    h(3,3) = 4;
    h(4,4) = 5;
    
    // Time interval
    double dt = 0.1;
    
    auto H = [h](vector<double> const& v){ return h*v; };
    
    // Initial vectors psi_0
    std::vector<vector<double>> psi0;
    
    // Eigenstate of H
    psi0.emplace_back(5);
    psi0[0](0) = 1;
    psi0[0](1) = 0.0;
    psi0[0](2) = 0.0;
    psi0[0](3) = 0.0;
    psi0[0](4) = 0.0;
    // Mixture of 2 eigenstates
    psi0.emplace_back(5);
    psi0[1](0) = 1.0/sqrt(3.0);
    psi0[1](1) = sqrt(2.0/3.0);
    psi0[1](2) = 0.0;
    psi0[1](3) = 0.0;
    psi0[1](4) = 0.0;
    // Mixture of 3 eigenstates
    psi0.emplace_back(5);
    psi0[2](0) = 1.0/sqrt(6.0);
    psi0[2](1) = 1.0/sqrt(3.0);;
    psi0[2](2) = 1.0/sqrt(2.0);;
    psi0[2](3) = 0.0;
    psi0[2](4) = 0.0;
    // Mixture of 4 eigenstates
    psi0.emplace_back(5);
    psi0[3](0) = 1.0/sqrt(10.0);
    psi0[3](1) = 1.0/sqrt(5.0);
    psi0[3](2) = sqrt(3.0/10.0);
    psi0[3](3) = sqrt(2.0/5.0);
    psi0[3](4) = 0.0;
    // Mixture of all 5 eigenstates
    psi0.emplace_back(5);
    psi0[4](0) = 1.0/sqrt(15.0);
    psi0[4](1) = sqrt(2.0/15.0);
    psi0[4](2) = 1.0/sqrt(5.0);
    psi0[4](3) = 2.0/sqrt(15.0);
    psi0[4](4) = 1.0/sqrt(3.0);
        
#ifdef KRYLOV_STATS
    krylov_params kp({10,1e-10,"krylov.stats.dat"});
#else
    krylov_params kp({10,1e-10});
#endif
    
    krylov_worker<decltype(H), vector<double>> kw(H,kp);
    
    std::vector<double> alpha, beta;
    for(int n = 0; n < 5; ++n){ 
        std::tie(alpha,beta) = kw(psi0[n]);
        // Check dimensions of Krylov's subspaces
        if(alpha.size() != n+1 || beta.size() != n) return -1;
        kw.reset();
    }

    // Final vectors psi (reference)
    std::vector<vector<double>> psi;

    psi.emplace_back(5);
    psi[0](0) = psi0[0](0)*exp(-h(0,0)*dt);
    psi[0](1) = 0.0;
    psi[0](2) = 0.0;
    psi[0](3) = 0.0;
    psi[0](4) = 0.0;

    psi.emplace_back(5);
    psi[1](0) = psi0[1](0)*exp(-h(0,0)*dt);
    psi[1](1) = psi0[1](1)*exp(-h(1,1)*dt);;
    psi[1](2) = 0.0;
    psi[1](3) = 0.0;
    psi[1](4) = 0.0;
    
    psi.emplace_back(5);
    psi[2](0) = psi0[2](0)*exp(-h(0,0)*dt);
    psi[2](1) = psi0[2](1)*exp(-h(1,1)*dt);
    psi[2](2) = psi0[2](2)*exp(-h(2,2)*dt);
    psi[2](3) = 0.0;
    psi[2](4) = 0.0;
    
    psi.emplace_back(5);
    psi[3](0) = psi0[3](0)*exp(-h(0,0)*dt);
    psi[3](1) = psi0[3](1)*exp(-h(1,1)*dt);
    psi[3](2) = psi0[3](2)*exp(-h(2,2)*dt);
    psi[3](3) = psi0[3](3)*exp(-h(3,3)*dt);
    psi[3](4) = 0.0;
    
    psi.emplace_back(5);
    psi[4](0) = psi0[4](0)*exp(-h(0,0)*dt);
    psi[4](1) = psi0[4](1)*exp(-h(1,1)*dt);
    psi[4](2) = psi0[4](2)*exp(-h(2,2)*dt);
    psi[4](3) = psi0[4](3)*exp(-h(3,3)*dt);
    psi[4](4) = psi0[4](4)*exp(-h(4,4)*dt);
    
    exp_h_worker<decltype(H), vector<double>> ehw(H,kp);

    for(int n = 0; n < 5; ++n){
        assert_all_close(ehw(psi0[n],dt), psi[n], 1e-10);
    }
    
    return 0;
}
