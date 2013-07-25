
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
#include <limits>
#include <complex>
#include <tuple>

#ifdef KRYLOV_STATS
#include "statistics.hpp"
#endif

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

struct krylov_params {
        
    std::size_t max_dim;
    double min_beta_threshold;
#ifdef KRYLOV_STATS
    std::string stats_file;
#endif
    
    static const std::size_t default_max_dim;
    static const double default_min_beta_threshold;
};

const std::size_t krylov_params::default_max_dim = std::numeric_limits<std::size_t>::max();
const double krylov_params::default_min_beta_threshold = 1e-7;
    
template <typename OperatorType, typename StateType>
    class krylov_worker {

    OperatorType const& H;

    typedef typename StateType::value_type scalar_type;

    std::vector<scalar_type> alpha;   // diagonal matrix elements
    std::vector<scalar_type> beta;    // subdiagonal matrix elements
   
    // Krylov basis states
    // H \approx V * T * V^+
    // Elements of 'basisstates' are columns of V
    std::vector<StateType> basisstates;
   
    // The tridiagonal matrix T has the following form
    // | alpha[0]    beta[0]     0       ...     |
    // | beta[0]     alpha[1]    beta[1]     ... |
    // | 0           beta[1]     alpha[2]    ... |
    // |             ...                         |
      
    // Temporaries
    StateType res_vector;
    scalar_type last_beta;
            
    static constexpr unsigned int reserved_krylov_dim = 5;
    
    // Adjustable parameters of the algorithm
    krylov_params kp;
    
#ifdef KRYLOV_STATS
    krylov_stats_collector stats;
#endif

#ifdef KRYLOV_STATS
    // Hash a pair of natural numbers (n,m)
    // n = 1, 2, ...
    // m = 1, 2, ..., n
    struct dims_hash {
        std::size_t operator()(std::pair<std::size_t,std::size_t> const& n_m) const
        {
            return n_m.first*(n_m.first-1)/2 + n_m.second - 1;
        }
    };
    
    std::unordered_map<std::pair<std::size_t,std::size_t>, std::size_t, dims_hash> dims_stats;
#endif
    
    // Calculates the next state in Krylov's basis.
    // Returns false if the previous state was an eigenstate of H
    void advance()
    {           
        basisstates.push_back(res_vector/last_beta);
        // try optimisation : if ok, use assign delegation to maintain genericty
	//res_vector.amplitudes()()=0; H.apply(basisstates.back(),res_vector);
        res_vector = H(basisstates.back());
        alpha.push_back(dotc(basisstates.back(),res_vector));
        res_vector -= alpha.back() * basisstates.back();
        if(beta.size() > 0) res_vector -= last_beta * basisstates[basisstates.size()-2];
        last_beta = std::sqrt(dotc(res_vector,res_vector));
        beta.push_back(last_beta);
    }
    
    public :

    typedef StateType state_type;
        
    krylov_worker(OperatorType const& H, krylov_params kp) :
        H(H), kp(kp)
#ifdef KRYLOV_STATS
        , stats(kp.stats_file)
#endif
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
    melements_t operator()(StateType const& initial_state)
    {        
        std::size_t space_dim = get_space_dim(initial_state);
        std::size_t max_dim = std::min(space_dim, kp.max_dim);
        reset();
        
        last_beta = 1.0;
        res_vector = initial_state;

        while(alpha.size() < max_dim && std::abs(last_beta) > kp.min_beta_threshold) advance();
        if(beta.size() != 0) beta.pop_back();

#ifdef KRYLOV_STATS
        stats(space_dim,alpha.size());
#endif
        
        return std::make_tuple(std::cref(alpha), std::cref(beta));
    }

    void reset()
    {
        alpha.clear();
        beta.clear();
        basisstates.clear();
    }

    template<typename KrylovCoeffs>
    StateType krylov_2_fock(KrylovCoeffs const& phi)
    {
        //StateType fock_state(basisstates[0].size());
        StateType st(res_vector); 
        //StateType st(res_vector.get_hilbert());
        //st() = 0;
        for(std::size_t i = 0; i < phi.size(); ++i)
            st += phi(i) * basisstates[i];  
        return st;
    }
    
#ifdef KRYLOV_STATS
    ~krylov_worker()
    {        
        stats.dump();
    }
#endif
    
  };

}}}}
#endif
