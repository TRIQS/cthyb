
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by M. Ferrero, O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef TRIQS_CTQMC_KRYLOV_KRYLOV_WORKER_HPP
#define TRIQS_CTQMC_KRYLOV_KRYLOV_WORKER_HPP
#include <vector>

#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/stev.hpp>

using triqs::arrays::blas::tridiag_worker;

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {
    
template <typename OperatorType, typename StateType>
    class krylov_worker {

    OperatorType const& H;

    typedef typename StateType::value_type scalar_type;
   
    // Krylov basis states
    // H \approx V * T * V^+
    // Elements of 'basisstates' are columns of V
    std::vector<StateType> basisstates;
   
    // The tridiagonal matrix T has the following form
    // | alpha[0]    beta[0]     0       ...     |
    // | beta[0]     alpha[1]    beta[1]     ... |
    // | 0           beta[1]     alpha[2]    ... |
    // |             ...                         |

    std::vector<scalar_type> alpha;   // diagonal matrix elements
    std::vector<scalar_type> beta;    // subdiagonal matrix elements
    
    // Temporaries
    StateType res_vector;
            
    static constexpr unsigned int reserved_krylov_dim = 5;
    
    // Adjustable parameters of the algorithm
    double gs_energy_convergence;
    
    // Tridiagonal matrix diagonalizer
    tridiag_worker tdw;
    
    // Returns the only matrix element of the 1x1 Krylov-projected matrix
    double first_iteration(StateType const& initial_state)
    {
        basisstates.push_back(initial_state);
        res_vector = H(basisstates.back());
        alpha.push_back(dotc(initial_state,res_vector));
        res_vector -= alpha.back() * initial_state;
        
        return alpha.back();
    }
    
    // Calculates the next state in Krylov's basis.
    // Returns false if the previous state was an eigenstate of H
    bool advance()
    {
        double new_beta = std::sqrt(dotc(res_vector,res_vector));
        // We don't really want to divide by zero
        if(std::abs(new_beta) < gs_energy_convergence) return false;
        
        beta.push_back(new_beta);
        basisstates.push_back(res_vector/new_beta);
        // try optimisation : if ok, use assign delegation to maintain genericty
        //res_vector.amplitudes()()=0; H.apply(basisstates.back(),res_vector);
        res_vector = H(basisstates.back());
        alpha.push_back(dotc(basisstates.back(),res_vector));
        res_vector -= alpha.back() * basisstates.back();
        res_vector -= beta.back() * basisstates[basisstates.size()-2];
    
        return true;
    }
    
    public :

    typedef StateType state_type;
        
    krylov_worker(OperatorType const& H, double gs_energy_convergence = 1e-10) :
        H(H), gs_energy_convergence(gs_energy_convergence),
        tdw(reserved_krylov_dim)
    {
        alpha.reserve(reserved_krylov_dim);
        beta.reserve(reserved_krylov_dim-1);
        basisstates.reserve(reserved_krylov_dim);
    };
    krylov_worker(krylov_worker const&) = default;
    krylov_worker& operator=(krylov_worker const&) = delete;

    // (main diagonal,subdiagonal) tuple
    typedef std::tuple< std::vector<scalar_type> const&,
                        std::vector<scalar_type> const&> melements_t;
    
    // initial_state MUST be of norm 1
    void operator()(StateType const& initial_state)
    {        
        reset();

        // First iteration
        double gs_energy = first_iteration(initial_state);
        tdw(alpha,beta);    // FIXME: no need to call this... but otherwise tdw.values() can not be called
        
        while(advance()){
            tdw(alpha,beta);
            if(std::abs(tdw.values()[0] - gs_energy) < gs_energy_convergence) break;
            gs_energy = tdw.values()[0];
        }
    }
    
    // Access eigenvalues and eigenvectors of the Krylov-projected operator
    arrays::vector_view<double> values() const { return tdw.values();}
    arrays::matrix_view<double> vectors() const { return tdw.vectors();}

    void reset()
    {
        alpha.clear();
        beta.clear();
        basisstates.clear();
    }

    template<typename KrylovCoeffs>
    StateType krylov_2_fock(KrylovCoeffs const& phi)
    {
        StateType st = make_zero_state(res_vector); 
        for(std::size_t i = 0; i < phi.size(); ++i)
            st += phi(i) * basisstates[i];
        return st;
    }
    
  };

}}}}
#endif
