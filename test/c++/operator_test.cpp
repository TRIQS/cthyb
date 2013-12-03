/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by I. Krivenko
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
#include <triqs/utility/first_include.hpp>

#include <iostream>
#include <vector>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <sstream>

#include "./operator.hpp"

using namespace triqs::utility;

int main(int argc, char **argv)
{
    // Operators without indices
    auto op_with_no_indices = c() + c_dag() - n();
    std::cout << "op_with_no_indices = " << op_with_no_indices << std::endl;
    
    // Operators with many indices
    auto op_with_many_indices = c(1,0.2,"a",true,-2) +
                            c_dag(3,0.15,"b",false,-5);
    std::cout << "op_with_many_indices = " << op_with_many_indices << std::endl;
                            
    // Commutation relations
#ifndef TRIQS_WORKAROUND_INTEL_COMPILER_BUGS
    std::vector<many_body_operator<double,int> > C = {c(1), c(2), c(3)};
    std::vector<many_body_operator<double,int> > Cd = {c_dag(1), c_dag(2), c_dag(3)};
#else
    std::vector<many_body_operator<double,int> > C; C.push_back(c(1)); C.push_back(c(3));C.push_back(c(3));
    std::vector<many_body_operator<double,int> > Cd; C.push_back(c_dag(1)); C.push_back(c_dag(3));C.push_back(c_dag(3));
#endif    
    std::cout << std::endl << "Anticommutators:" << std::endl;
    for(auto const& cdi : Cd)
    for(auto const& ci : C){
        std::cout << "{" << cdi << ", " << ci << "} = " << cdi*ci + ci*cdi << std::endl;
    }
    
    std::cout << std::endl << "Commutators:" << std::endl;
    for(auto const& cdi : Cd)
    for(auto const& ci : C){
        std::cout << "[" << cdi << ", " << ci << "] = " << cdi*ci - ci*cdi << std::endl;
    }
    
    // Algebra
    auto x = c(0);
    auto y = c_dag(1);

    std::cout << std::endl << "Algebra:" << std::endl;    
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;

    std::cout << "-x = " << -x << std::endl;
    std::cout << "x + 2.0 = " << x + 2.0 << std::endl;
    std::cout << "2.0 + x = " << 2.0 + x << std::endl;
    std::cout << "x - 2.0 = " << x - 2.0 << std::endl;
    std::cout << "2.0 - x = " << 2.0 - x << std::endl;
    std::cout << "3.0*y = " << 3.0*y << std::endl;
    std::cout << "y*3.0 = " << y*3.0 << std::endl;
    std::cout << "x + y = " << x + y << std::endl;
    std::cout << "x - y = " << x - y << std::endl;
    std::cout << "(x + y)*(x - y) = " << (x + y)*(x - y) << std::endl;

    // N^3
    std::cout << std::endl << "N^3:" << std::endl;
    auto N = n("up") + n("dn");
    auto N3 = N*N*N;
    std::cout << "N = " << N << std::endl;
    std::cout << "N^3 = " << N3 << std::endl; 
        
    // Serialization
    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa & N3;

    boost::archive::text_iarchive ia(ss);
    many_body_operator<double,std::string> new_N3;
    ia & new_N3;
    
    std::cout << "New N^3 = " << new_N3 << std::endl; 
    
    return 0;
}
