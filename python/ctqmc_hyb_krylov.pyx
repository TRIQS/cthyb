from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from pytriqs.gf.local.gf cimport *
from pytriqs.parameters.parameters cimport *
from many_body_operator cimport *

import numpy
import pytriqs.utility.mpi as mpi
from pytriqs.gf.local import *

ctypedef many_body_operator[double,string,string] operator_c

cdef extern from "c++/fundamental_operator_set.hpp" namespace "triqs::app::impurity_solvers::ctqmc_krylov":
    cdef cppclass fundamental_operator_set[string,string]:
        
        fundamental_operator_set() except +
        void add_operator(string, string)
        
cdef extern from "c++/sorted_spaces.hpp" namespace "triqs::app::impurity_solvers::ctqmc_krylov":
    cdef cppclass block_desc_t[string,string]:
        
        string name
        vector[pair[string,string]] indices
        
        block_desc_t() except +

cdef extern from "c++/ctqmc_krylov.hpp" namespace "triqs::app::impurity_solvers":

    cdef cppclass solver_c "triqs::app::impurity_solvers::ctqmc_krylov::ctqmc_krylov":
        
      solver_c(parameters p,
              const operator_c& h_loc,
              const vector[operator_c] & quantum_numbers,
              fundamental_operator_set & fops,
              const vector[block_desc_t] & block_stucture
              ) except +

      #container views

      #input containers
      gf_block_imtime & deltat_view()
      #gf_block_imfreq & g0w_view()

      #imaginary-time Green's functions
      gf_block_imtime & gt_view()

      #Matsubara Green's functions
      #gf_block_imfreq & gw_view()

      void solve(parameters p) except +

      parameter_defaults constructor_defaults() 
      parameter_defaults solve_defaults() 
  
cdef extern from "triqs/utility/formatted_output.hpp" namespace "triqs::utility":
  
    cdef string print_formatted( vector[vector[std_string]] &) except +

cdef class Solver:

    cdef solver_c * _c
    cdef object block_indices_pack

    def __cinit__(self, **kw):
        
        # Parameters
        p = kw['parameters']
        
        # Hamiltonian
        cdef operator_c H_local = (<Operator?> kw['H_local'])._c
        
        # Quantum numbers
        cdef vector[operator_c] quantum_numbers
        for qn in kw['quantum_numbers']: quantum_numbers.push_back((<Operator?> qn)._c)
        
        cdef vector[block_desc_t] block_stucture
        cdef fundamental_operator_set fops
        
        self.block_indices_pack = []
        
        cdef string block_name, index_name
        cdef block_desc_t block
        for b_name in kw['gf_struct']:    # Legacy name
            block_indices = kw['gf_struct'][b_name]
            block_name = str(b_name)
            block.name = block_name
            block.indices.clear()
            
            for i_name in block_indices:
                index_name = str(i_name)    
                fops.add_operator(block_name,index_name)
                block.indices.push_back(pair[string,string](block_name,index_name))
            
            block_stucture.push_back(block)
            self.block_indices_pack.append([range(block.indices.size()),range(block.indices.size())])
        
        self._c = new solver_c((<Parameters?> p)._c,H_local,quantum_numbers,fops,block_stucture)

    def __dealloc__(self):
        del self._c

    def solve(self, **kw):
        self._c.solve((<Parameters?>kw['parameters'])._c)

    #accessors for input containers (read-write)
    property Delta_tau:
        """Hybridization function"""
        def __get__(self): return make_BlockGfImTime(self._c.deltat_view(), self.block_indices_pack, "Delta(tau)")
        def __set__(self,val):
            cdef g0 = make_BlockGfImTime(self._c.deltat_view(), self.block_indices_pack, "Delta(tau)")
            g0 <<= val

    #accessors for output containers (read-only)    
    property G_tau:
        """Imaginary-time Green's function"""
        def __get__(self): return make_BlockGfImTime(self._c.gt_view(), self.block_indices_pack, "G") # backward compatibility to h5
        #def __get__(self): return make_BlockGfImTime(self._c.gt_view(), self.block_indices_pack, "G(tau)")

    # Todo : add option to generate correct rst table..
    def help(self) :
        """Generate the documentation of the solver"""
        #cdef vector[vector[std_string]] h
        s = """ 
Parameters of the Krylov CT-HYB solver :
        
Constructor :
%s 
        
Solve method : 
%s
"""% (print_formatted (self._c.constructor_defaults().generate_help())
         ,print_formatted (self._c.solve_defaults().generate_help())
         )

        return s