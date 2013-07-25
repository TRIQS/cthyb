from libcpp.string cimport string
from libcpp cimport bool
from arrays cimport *
from pytriqs.gf.local.gf cimport *
import numpy

from pytriqs.gf.local import *
import pytriqs.utility.mpi as mpi
from pytriqs.parameters.parameters cimport *

cdef extern from "c++/ctqmc_seg.hpp" namespace "triqs::app::impurity_solvers":

    cdef cppclass solver_c "triqs::applications::impurity_solvers::ctqmc_seg::ctqmc_seg":
      solver_c(parameters &) except +

      #container views

      #input containers
      gf_block_imtime & deltat_view()
      gf_block_imfreq & g0w_view()
      gf_imtime & kt_view()
      gf_imtime & kprimet_view()

      #imaginary-time Green's functions
      gf_block_imtime & gt_view()
      gf_block_imtime & ft_view()
      gf_imtime & nnt_view()

      #Matsubara Green's functions
      gf_block_imfreq & gw_view()
      gf_block_imfreq & fw_view()
      gf_block_imfreq & sw_view()
      gf_imfreq & nnw_view()

      #imaginary-time Green's functions
      gf_block_legendre & gl_view()
      gf_block_legendre & fl_view()

      #miscellaneous
      matrix_view[double] histogram()
      matrix_view[double] nn()

      void solve(parameters &) except +
      double percent_done()

      parameter_defaults constructor_defaults() 
      parameter_defaults solve_defaults() 
  
cdef extern from "triqs/utility/formatted_output.hpp" namespace "triqs::utility":
  
    cdef std_string print_formatted( vector[vector[std_string]] &) except +

cdef class Solver:

    cdef solver_c * _c
    cdef int Verbosity
    cdef object param    
    cdef object block_indices_pack
    cdef object block_indices_pack_bosonic
    cdef bool First

    def __init__(self, **kw):
        assert ('parameters' in kw and len(kw)==1) or ('parameters' not in kw), 'wrong arguments'
        self.param=kw['parameters'] if 'parameters' in kw else Parameters().update(kw)
        block_names=self.param['block_names']
        self.block_indices_pack = []
        for i in block_names:  self.block_indices_pack.append( [[i], [i]] )
        self.block_indices_pack_bosonic = [block_names, block_names] 
        self.Verbosity = 2 if mpi.rank ==0 else 0
        self._c = new solver_c((<Parameters?>self.param)._c)
        self.First = True    

    def __dealloc__(self):
        del self._c

    #accessors for input containers (read-write)
    property Delta_tau:
        """Hybridization function"""
        def __get__(self): return make_BlockGfImTime(self._c.deltat_view(), self.block_indices_pack, "Delta(tau)")
        def __set__(self,val):
            cdef g0 = make_BlockGfImTime(self._c.deltat_view(), self.block_indices_pack, "Delta(tau)")
            g0 <<= val

    property G0:
        """Bare Green's function"""
        def __get__(self): return make_BlockGfImFreq(self._c.g0w_view(), self.block_indices_pack, "G0(omega)")
        def __set__(self,val):
            cdef g0 = make_BlockGfImFreq(self._c.g0w_view(), self.block_indices_pack, "G0(omega)")
            g0 <<= val

    property K_tau:
        """Retarded interaction kernel"""
        def __get__(self): return make_GfImTime(self._c.kt_view(), self.block_indices_pack_bosonic, "K(tau)")
        def __set__(self,val):
            cdef g0 = make_GfImTime(self._c.kt_view(), self.block_indices_pack_bosonic, "K(tau)")
            g0 <<= val

    property Kprime_tau:
        """Derivative of retarded interaction kernel"""
        def __get__(self): return make_GfImTime(self._c.kprimet_view(), self.block_indices_pack_bosonic, "Kprime(tau)")
        def __set__(self,val):
            cdef g0 = make_GfImTime(self._c.kprimet_view(), self.block_indices_pack_bosonic, "Kprime(tau)")
            g0 <<= val

    #accessors for output containers (read-only)
    property nn_tau:
        """Density-density correlation function in imaginary time"""
        def __get__(self): return make_GfImTime(self._c.nnt_view(), self.block_indices_pack_bosonic, "nn(tau)")

    property nn_omega:
        """Density-density correlation function on Matsubara frequencies"""
        def __get__(self): return make_GfImFreq(self._c.nnw_view(), self.block_indices_pack_bosonic, "nn(omega)")
    
    property G_tau:
        """Imaginary-time Green's function"""
        def __get__(self): return make_BlockGfImTime(self._c.gt_view(), self.block_indices_pack, "G") # backward compatibility to h5
        #def __get__(self): return make_BlockGfImTime(self._c.gt_view(), self.block_indices_pack, "G(tau)")

    property F_tau:
        """Improved estimator in imaginary time"""
        def __get__(self): return make_BlockGfImTime(self._c.ft_view(), self.block_indices_pack, "F(tau)")

    property G_legendre:
        """Imaginary-time Green's function"""
        def __get__(self): return make_BlockGfLegendre(self._c.gl_view(), self.block_indices_pack, "G(l)")

    property F_legendre:
        """Improved estimator in imaginary time"""
        def __get__(self): return make_BlockGfLegendre(self._c.fl_view(), self.block_indices_pack, "F(l)")

    property G_omega:
        """Matsubara Green's function"""
        def __get__(self): return make_BlockGfImFreq(self._c.gw_view(), self.block_indices_pack, "G(omega)")

    property F_omega:
        """Improved estimator on Matsubara frequencies"""
        def __get__(self): return make_BlockGfImFreq(self._c.fw_view(), self.block_indices_pack, "F(omega)")

    property Sigma_omega:
        """Self-energy on Matsubara frequencies"""
        def __get__(self): return make_BlockGfImFreq(self._c.sw_view(), self.block_indices_pack, "Sigma(omega)")

    property hist:
        """Histogram of pert order"""
        def __get__(self): return self._c.histogram().to_python()

    property nn:
        """Density-density static correlations"""
        def __get__(self): return self._c.nn().to_python()

    def solve(self, **kw):
        """ Solve the impurity problem """
        assert ('parameters' in kw and len(kw)==1) or ('parameters' not in kw), 'wrong arguments'
        self.param=Parameters().update(kw)
        #self.param=self.param.update(kw)
        # Call the solver
        self._c.solve((<Parameters?>self.param)._c)

    def percent_done(self): return self._c.percent_done()

    # Todo : add option to generate correct rst table..
    def help(self) :
        """Generate the documentation of the solver"""
        #cdef vector[vector[std_string]] h
        s = """ 
Parameters of the cthyb segment solver :
        
Constructor :
%s 
        
Solve method : 
%s
"""% (print_formatted (self._c.constructor_defaults().generate_help())
         ,print_formatted (self._c.solve_defaults().generate_help())
         )

        return s
