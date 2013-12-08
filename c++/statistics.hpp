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
#ifndef TRIQS_CTQMC_STATISTICS_HPP
#define TRIQS_CTQMC_STATISTICS_HPP

#include <unordered_map>
#include <fstream>
#include <tuple>
#include <boost/mpi.hpp>
#include <boost/serialization/hash_collections_load_imp.hpp>
#include <boost/serialization/hash_collections_save_imp.hpp>

namespace cthyb_krylov {

    
struct dims_stats_collector {
    
    // Hash a pair of natural numbers (n,m)
    // n = 1, 2, ...
    // m = 1, 2, ..., n
    struct dims_hash {
        std::size_t operator()(std::pair<std::size_t,std::size_t> const& n_m) const
        {
            return n_m.first*(n_m.first-1)/2 + n_m.second - 1;
        }
    };
    
    std::string stats_file;
    typedef std::unordered_map<std::pair<std::size_t,std::size_t>, unsigned long, dims_hash> dims_stats_t;
    dims_stats_t dims_stats;
    
public:
    
    explicit dims_stats_collector(std::string const& stats_file) : dims_stats(), stats_file(stats_file) {}
    
    void operator()(std::size_t n, std::size_t m, unsigned long count = 1)
    {
        bool inserted;
        dims_stats_t::iterator it;
        std::tie(it,inserted) = dims_stats.insert(std::make_pair(std::make_pair(n,m),count));
        if(!inserted) it->second += count;
    }
    
    void dump()
    {
        bool write_file = true;
        if(boost::mpi::environment::initialized()){
            boost::mpi::communicator c;
            if(c.rank() == 0){
                std::size_t counts_to_recv;
                boost::mpi::reduce(c, dims_stats.size(), counts_to_recv, std::plus<std::size_t>(), 0);
                counts_to_recv -= dims_stats.size();
                
                for(;counts_to_recv>0; --counts_to_recv){
                    dims_stats_t::value_type v;
                    c.recv(boost::mpi::any_source, 10000, v);
                    (*this)(v.first.first,v.first.second,v.second);
                }

            }else{
                boost::mpi::reduce(c, dims_stats.size(), std::plus<std::size_t>(), 0);
                for(auto const& v : dims_stats) c.send(0, 10000, v);
                write_file = false;
            }
        }
        
        if(write_file){
            std::ofstream file(stats_file);
            file << "# space_dim\teffective_dim\tcount" << std::endl;
            for(auto const& s : dims_stats)
                file << s.first.first << '\t' << s.first.second << '\t' << s.second << std::endl;
            file.close();
        }
    }
    
};

}

#endif
