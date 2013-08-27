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

template<typename HamiltonianType, typename StateType, bool UseKrylov = true>
class exp_h_worker {};

// Use Krylov subspace projection to calculate \exp(-\tau*H)
template<typename HamiltonianType, typename StateType>
class exp_h_worker<HamiltonianType, StateType, true> {
    
    krylov_worker<HamiltonianType,StateType> kw;
    
public:
    
    exp_h_worker(HamiltonianType const& H, krylov_params kp) : kw(H,kp) {};
    exp_h_worker(exp_h_worker const&) = default;
    exp_h_worker& operator=(exp_h_worker const&) = delete;
    
    typedef StateType state_type;
    typedef typename state_type::value_type scalar_type;
    
    state_type operator()(state_type const& initial_state, double dtau)
    {
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
    }
};

// Use direct matrix-matrix multiplication to calculate \exp(-\tau*H)
template<typename HamiltonianType, typename StateType>
class exp_h_worker<HamiltonianType, StateType, false> {
    
public:
    
    exp_h_worker() {};
    exp_h_worker(exp_h_worker const&) = default;
    exp_h_worker& operator=(exp_h_worker const&) = delete;
    
    typedef StateType state_type;
    typedef typename state_type::value_type scalar_type;
        
    state_type operator()(state_type const& initial_state, double dtau, sorted_spaces::eigensystem_t const& eigensystem)
    {
        auto const& eigenvalues = eigensystem.eigenvalues;
        auto const& unitary_matrix = eigensystem.unitary_matrix;
        std::size_t dim = eigenvalues.size();
        
        matrix<scalar_type> matrix_exp(dim,dim);
        matrix_exp() = 0;
        for(std::size_t n = 0; n < dim; ++n)
            matrix_exp(n,n) = exp(-dtau*(eigensystem.eigenvalues(n)));
        
        matrix_exp = unitary_matrix * matrix_exp * unitary_matrix.transpose();

        StateType st = make_zero_state(initial_state);
        // FIXME: not supposed to work with the map-based version of state...
        st.amplitudes() = matrix_exp * initial_state.amplitudes();
        
        return st;
    }
};
    
}}}}
#endif
