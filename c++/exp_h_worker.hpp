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

#include <tuple>
#include <algorithm>

#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/stev.hpp>
#include "krylov_worker.hpp"

using namespace triqs::arrays;
using triqs::arrays::blas::tridiag_worker;

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

template<typename HamiltonianType, typename StateType>
class exp_h_worker {
    
    typedef krylov_worker<HamiltonianType,StateType> kw_type;
    kw_type kw;
    
public:
    
    exp_h_worker(HamiltonianType const& H, krylov_params kp) : kw(H,kp) {};
    exp_h_worker(exp_h_worker const&) = default;
    exp_h_worker& operator=(exp_h_worker const&) = delete;
    
    typedef StateType state_type;
    typedef typename state_type::value_type scalar_type;
    
    state_type operator()(state_type const& initial_state, double dtau)
    {
        scalar_type initial_state_norm  = std::sqrt(dotc(initial_state,initial_state));
        auto melements = kw(initial_state/initial_state_norm);
    
        std::size_t krylov_dim = std::get<0>(melements).size();
        auto tdw = tridiag_worker(2*krylov_dim);
    
        // Diagonalize
        tdw(std::get<0>(melements),std::get<1>(melements));
        auto eigenvalues = tdw.values();
        
        matrix<scalar_type> krylov_exp(krylov_dim,krylov_dim);
        krylov_exp() = 0;
        for(std::size_t n = 0; n < krylov_dim; ++n)
            krylov_exp(n,n) = exp(-dtau*eigenvalues(n));
    
        // FIXME: Hermitian adjoint instead of transpose?
        krylov_exp = tdw.vectors().transpose() * krylov_exp * tdw.vectors();
            
        auto krylov_coeffs = initial_state_norm * krylov_exp(ellipsis(),0);
        
        return kw.krylov_2_fock(krylov_coeffs);
    }
};
    
}}}}
#endif
