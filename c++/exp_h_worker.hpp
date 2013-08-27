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
#ifndef TRIQS_CTQMC_KRYLOV_EXP_H_HPP
#define TRIQS_CTQMC_KRYLOV_EXP_H_HPP

#include "krylov_worker.hpp"
#include "sorted_spaces.hpp"

using namespace triqs::arrays;

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

template<typename HamiltonianType, typename StateType>
class exp_h_worker {
    
    krylov_worker<HamiltonianType,StateType> kw;
    sorted_spaces sosp;
    
    std::size_t small_matrix_size;
    
public:
    
    exp_h_worker(HamiltonianType const& H, sorted_spaces const& sosp, krylov_params kp, std::size_t small_matrix_size) :
        kw(H,kp), sosp(sosp), small_matrix_size(small_matrix_size) {};
    exp_h_worker(exp_h_worker const&) = default;
    exp_h_worker& operator=(exp_h_worker const&) = delete;
    
    typedef StateType state_type;
    typedef typename state_type::value_type scalar_type;
    
    state_type operator()(state_type const& initial_state, double dtau)
    {
        auto const& space = initial_state.get_hilbert();
        std::size_t space_dim = space.dimension();
        
        if(space_dim > small_matrix_size){
        
            scalar_type initial_state_norm  = std::sqrt(dotc(initial_state,initial_state));
            kw(initial_state/initial_state_norm);

            auto eigenvalues = kw.values();
            std::size_t krylov_dim = eigenvalues.size();
        
            matrix<scalar_type> krylov_exp(krylov_dim,krylov_dim);
            krylov_exp() = 0;
            for(std::size_t n = 0; n < krylov_dim; ++n)
                krylov_exp(n,n) = exp(-dtau*eigenvalues(n));
        
            krylov_exp = kw.vectors().transpose() * krylov_exp * kw.vectors();
            auto krylov_coeffs = initial_state_norm * krylov_exp(ellipsis(),0);
   
            return kw.krylov_2_fock(krylov_coeffs);
            
        } else {
            
            auto const& eigensystem = sosp.get_eigensystems()[space.get_index()];
            auto const& eigenvalues = eigensystem.eigenvalues;
            auto const& unitary_matrix = eigensystem.unitary_matrix;
        
            matrix<scalar_type> matrix_exp(space_dim,space_dim);
            matrix_exp() = 0;
            for(std::size_t n = 0; n < space_dim; ++n)
                matrix_exp(n,n) = exp(-dtau*(eigensystem.eigenvalues(n)));
        
            matrix_exp = unitary_matrix * matrix_exp * unitary_matrix.transpose();

            StateType st = make_zero_state(initial_state);
            // FIXME: not supposed to work with the map-based version of state...
            st.amplitudes() = matrix_exp * initial_state.amplitudes();
        
            return st;
        }
    }
};
    
}}}}
#endif
