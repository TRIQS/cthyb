#!/bin/env pytriqs

import time
from pytriqs.gf.local import *
from pytriqs.operators import *
from pytriqs.operators.util.op_struct import set_operator_structure, get_mkind
from pytriqs.operators.util.U_matrix import U_matrix
from pytriqs.operators.util.hamiltonians import h_int_kanamori, h_int_slater
from pytriqs.applications.impurity_solvers.cthyb import *
import numpy as np

spin_names = ("up","dn")

def partition(h_int,h_k,gf_struct,QN=None):
    p = {}
    p["verbosity"] = 0
    p["max_time"] = -1
    p["length_cycle"] = 1
    p["n_warmup_cycles"] = 0
    p["n_cycles"] = 0

    if QN is None:
        p['partition_method'] = "autopartition"
    else:
        p['partition_method'] = "quantum_numbers"
        p['quantum_numbers'] = QN

    S = SolverCore(beta=1, gf_struct=gf_struct)
    S.G0_iw << inverse(iOmega_n - h_k)

    start = time.clock()
    S.solve(h_int=h_int,**p)
    end = time.clock()

    return S.eigensystems, end-start

table_format = ("%40s "+"%20s "*3)
print table_format % ("Model","Dimension of HS","Quantum numbers","Autopartition")
print table_format % ("=====","===============","===============","=============")

def print_line(model,hs_size,qn,ap):
    print table_format % (model,hs_size,"%i (%.3f sec)" % qn,"%i (%.3f sec)" % ap)

### Kanamori Hamiltonians
def run_kanamori(max_orbitals,orbital_mixing):
    for num_orbitals in range(2,max_orbitals+1):
        orb_names = range(num_orbitals)
        gf_struct = set_operator_structure(spin_names,orb_names,True)
        mkind = get_mkind(True,None)

        U = 1.0*(np.ones((num_orbitals,num_orbitals))-np.eye(num_orbitals))
        Up = 2.0*np.ones((num_orbitals,num_orbitals))
        J = 0.2
        V = 0.3

        h_k = V*np.ones((num_orbitals,num_orbitals)) if orbital_mixing else V*np.eye(num_orbitals)
        h_int = h_int_kanamori(spin_names,orb_names,U,Up,J,True)

        # Quantum numbers
        QN = [sum([n(*mkind("up",o)) for o in orb_names],Operator()),   # N_up
              sum([n(*mkind("dn",o)) for o in orb_names],Operator())]   # N_down
        if not orbital_mixing:
            # PS quantum number
            QN.append(Operator())
            for i, o in enumerate(orb_names):
                dn = n(*mkind("up",o)) - n(*mkind("dn",o))
                QN[2] += (2**i)*dn*dn

        eig_qn,time_qn = partition(h_int,h_k,gf_struct,QN)
        eig_ap,time_ap = partition(h_int,h_k,gf_struct)

        model = "Kanamori, %i orbitals"%num_orbitals
        if orbital_mixing: model += " (orbital mixing)"
        print_line(model,2**(2*num_orbitals),(len(eig_qn),time_qn),(len(eig_ap),time_ap))

### Slater Hamiltonians
def run_slater(L,is_cubic):
    for l in L:
        orb_names = range(2*l+1)
        num_orbitals = len(orb_names)
        gf_struct = set_operator_structure(spin_names,orb_names,True)
        mkind = get_mkind(True,None)

        F0, F2, F4 = 3.0, 0.6, 0.1

        U_mat = U_matrix(l,[F0,F2,F4],basis='cubic' if is_cubic else 'spherical')
        h_int = h_int_slater(spin_names,orb_names,U_mat,True)

        h_k = np.zeros((num_orbitals,num_orbitals))

        # Quantum numbers
        QN = [sum([n(*mkind("up",o)) for o in orb_names],Operator()),   # N_up
              sum([n(*mkind("dn",o)) for o in orb_names],Operator())]   # N_down

        eig_qn,time_qn = partition(h_int,h_k,gf_struct,QN)
        eig_ap,time_ap = partition(h_int,h_k,gf_struct)

        model = "Slater, %i orbitals"%num_orbitals
        model += (" (cubic basis)" if is_cubic else " (spherical basis)")
        print_line(model,2**(2*num_orbitals),(len(eig_qn),time_qn),(len(eig_ap),time_ap))

print
run_kanamori(7,False)
print
run_kanamori(7,True)
print
run_slater([2,3],False)
print
run_slater([2,3],True)