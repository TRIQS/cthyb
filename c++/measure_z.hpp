
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
#ifndef TRIQS_CTQMC_KRYLOV_MEASURE_Z_HPP
#define TRIQS_CTQMC_KRYLOV_MEASURE_Z_HPP

#include <triqs/utility/first_include.hpp>

#include <complex>
#include <functional>
#include <boost/mpi/collectives.hpp>

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {

// Measure for the partition function and average sign
struct measure_z {

    typedef std::complex<double> mc_sign_type;
    
    // the partition function and total number of measures
    mc_sign_type z;
    long long num;
    mc_sign_type average_sign;

    measure_z(): z(0), num(0) {}

    void accumulate(mc_sign_type s) { z += s; num += 1; }

    void collect_results(boost::mpi::communicator const &c) {
        boost::mpi::all_reduce(c, z, z, std::plus<mc_sign_type>());
        boost::mpi::all_reduce(c, num, num, std::plus<long long>());
        average_sign = z / num;
    }
  
};
  
}}}}

#endif 
