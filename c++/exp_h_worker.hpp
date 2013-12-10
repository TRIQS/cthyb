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

//#define KRYLOV_STATS

#ifdef KRYLOV_STATS
#include "statistics.hpp"
#define DIMS_STATS_FILE "krylov_stat.txt"
#define BOUNDARY_STATS_FILE "bs_stat.txt"
#endif

using namespace triqs::arrays;

namespace cthyb_krylov {

template<typename HamiltonianType, typename StateType>
class exp_h_worker {
    
    typedef StateType state_type;
    typedef typename state_type::scalar_t scalar_type;
    
    krylov_worker<HamiltonianType,StateType> kw;
    sorted_spaces sosp;
    
    int small_matrix_size;
    
    // Temporary matrices 
    matrix<scalar_type> matrix_exp;
    matrix<scalar_type> krylov_exp;
	triqs::arrays::vector<scalar_type> v2;


#ifdef KRYLOV_STATS
    dims_stats_collector stats;
#endif
    
public:
    
    exp_h_worker(HamiltonianType const& H, sorted_spaces const& sosp, double gs_energy_convergence, int small_matrix_size) :
        kw(H,gs_energy_convergence),
        sosp(sosp),
        small_matrix_size(small_matrix_size)
#ifdef KRYLOV_STATS
        ,stats(DIMS_STATS_FILE)
#endif
    {
        int  max_subspace_dim = 0;
        for(int nsp = 0; nsp < sosp.n_subspaces(); ++nsp)
            max_subspace_dim = std::max(max_subspace_dim,sosp.subspace(nsp).dimension());

        small_matrix_size = std::min(small_matrix_size,max_subspace_dim);
        
        matrix_exp.resize(small_matrix_size,small_matrix_size);
        krylov_exp.resize(max_subspace_dim,max_subspace_dim);
	v2.resize(max_subspace_dim);
    }

    exp_h_worker(exp_h_worker const&) = default;
    exp_h_worker& operator=(exp_h_worker const&) = delete;
        
    state_type operator()(state_type const& initial_state, double dtau)
    {
        auto const& space = initial_state.get_hilbert();
        std::size_t space_dim = space.dimension();
      
             if(space_dim > small_matrix_size){
        
            scalar_type initial_state_norm  = std::sqrt(dot_product(initial_state,initial_state));
            kw(initial_state/initial_state_norm);

            auto eigenvalues = kw.values();
            std::size_t krylov_dim = eigenvalues.size();

#ifdef KRYLOV_STATS
            stats(space_dim,krylov_dim);
#endif
            auto all = range(0,krylov_dim);
            
            krylov_exp(all,all) = 0;
            for(std::size_t n = 0; n < krylov_dim; ++n)
                krylov_exp(n,n) = exp(-dtau*eigenvalues(n));
        
            krylov_exp(all,all) = kw.vectors().transpose() * krylov_exp(all,all) * kw.vectors();
            auto krylov_coeffs = initial_state_norm * krylov_exp(all,0);
   
            return kw.krylov_2_fock(krylov_coeffs);
            
        } else {

#ifdef KRYLOV_STATS
            stats(space_dim,space_dim);
#endif
            
            auto const& eigensystem = sosp.get_eigensystems()[space.get_index()];
            auto const& eigenvalues = eigensystem.eigenvalues;
            auto const& unitary_matrix = eigensystem.unitary_matrix;
        
            auto all = range(0,space_dim);
            StateType st = make_zero_state(initial_state);
           
//#define NOT_OPTIMAL 
#ifdef NOT_OPTIMAL 
            matrix_exp(all,all) = 0;
            for(std::size_t n = 0; n < space_dim; ++n)
                matrix_exp(n,n) = exp(-dtau*(eigensystem.eigenvalues(n)));
        
            matrix_exp(all,all) = unitary_matrix * matrix_exp(all,all) * unitary_matrix.transpose();

            // FIXME: not supposed to work with the map-based version of state...
            
	    st.amplitudes() = matrix_exp(all,all) * initial_state.amplitudes();
#else
	    //v2(all) = unitary_matrix.transpose() * initial_state.amplitudes();
	    triqs::arrays::blas::gemv(1.0, unitary_matrix.transpose(),  initial_state.amplitudes(), 0.0, v2(all)) ;
            for(int n = 0; n < space_dim; ++n) v2[n] *= exp(-dtau* eigensystem.eigenvalues(n));
	    //std::cout  << "apply exp : min energy :"<<  eigensystem.eigenvalues(0) << "  tau = "<< dtau << std::endl; 
            //st.amplitudes() =  unitary_matrix * v2(all);		
	    triqs::arrays::blas::gemv(1.0, unitary_matrix, v2(all), 0.0,  st.amplitudes()) ;
#endif

            return st;
        }
    }
        void apply(state_type & initial_state, double dtau)
    {
        auto const& space = initial_state.get_hilbert();
        std::size_t space_dim = space.dimension();
      
            auto const& eigensystem = sosp.get_eigensystems()[space.get_index()];
            auto const& eigenvalues = eigensystem.eigenvalues;
            auto const& unitary_matrix = eigensystem.unitary_matrix;
        
            auto all = range(0,space_dim);
            //StateType st = make_zero_state(initial_state);
           
	    //v2(all) = unitary_matrix.transpose() * initial_state.amplitudes();
	    triqs::arrays::blas::gemv(1.0, unitary_matrix.transpose(),  initial_state.amplitudes(), 0.0, v2(all)) ;
            for(int n = 0; n < space_dim; ++n) v2[n] *= exp(-dtau* eigensystem.eigenvalues(n));
            //st.amplitudes() =  unitary_matrix * v2(all);		
	    triqs::arrays::blas::gemv(1.0, unitary_matrix, v2(all), 0.0,  initial_state.amplitudes()) ;

            //return st;
    }

#ifdef KRYLOV_STATS
    ~exp_h_worker()
    {        
        stats.dump();
    }
#endif

};
    
}
#endif
