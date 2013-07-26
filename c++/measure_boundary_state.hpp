
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by I. Krivenko, M. Ferrero, O. Parcollet
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
#ifndef TRIQS_CTQMC_KRYLOV_MEASURE_BOUNDRY_STATE_HPP
#define TRIQS_CTQMC_KRYLOV_MEASURE_BOUNDRY_STATE_HPP

#include <triqs/utility/first_include.hpp>

#include <vector>
#include <functional>
#include <algorithm>
#include <fstream>
#include <boost/serialization/vector.hpp>
#include <boost/mpi/collectives.hpp>

#include "qmc_data.hpp"

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

// Measure imaginary time Green's function (one block)
struct measure_boundary_state {
    
    typedef std::complex<double> mc_sign_type;
        
    qmc_data const& data;
    std::string stats_file;
    
    typedef std::vector<std::size_t> hist_t;
    std::vector<hist_t> stats;
    
    measure_boundary_state(qmc_data const& data, std::string const& stats_file) :
        data(data), stats_file(stats_file), stats(data.sosp.n_subspaces())
    {
        for(std::size_t spn=0; spn<stats.size(); ++spn){
            stats[spn].resize(data.sosp.subspace(spn).dimension(),0);
        }
    }

    void accumulate(mc_sign_type s) {
        // Accumulates statistics for the first element of boundary_block_states_ids only.
        for(std::size_t spn=0; spn<stats.size(); ++spn){
            stats[spn][data.config.boundary_block_states_ids[spn][0]]++;
        }
    }

    void collect_results(boost::mpi::communicator const &c) {
        
        auto sum_vector = [](hist_t const& x, hist_t const& y) 
            {
                hist_t result(x.size(),0);
                transform(x.cbegin(),x.cend(),y.cbegin(),result.begin(),std::plus<std::size_t>());
                return result;
            };
            
        for(std::size_t spn=0; spn<stats.size(); ++spn){
            boost::mpi::all_reduce(c,stats[spn],stats[spn],sum_vector);
        }
        
        if(c.rank()==0){
            std::ofstream file(stats_file);
            file << "# Eigensystems" << std::endl;
            for(std::size_t spn=0; spn<data.sosp.n_subspaces(); ++spn){
                auto const& es = data.sosp.get_eigensystems()[spn];
                file << "Subspace " << spn << std::endl;
                for(std::size_t n=0; n<es.eigenvalues.size(); ++n){
                    file << n << "\t" << es.eigenvalues[n] << "\t" << es.eigenstates[n] << std::endl;
                }
            }
            file << std::endl << "# Statistics" << std::endl;
            for(std::size_t spn=0; spn<data.config.boundary_block_states_ids.size(); ++spn){
                auto const& s = stats[spn];
                file << "Subspace " << spn << std::endl;
                for(std::size_t n=0; n<s.size(); ++n){
                    file << n << "\t" << s[n] << std::endl;
                }
            }            
        }
    }
  
};
  
}}}}

#endif
