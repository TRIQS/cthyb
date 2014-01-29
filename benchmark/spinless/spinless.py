#!/bin/env pytriqs

import numpy

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.parameters.parameters import Parameters
from pytriqs.applications.impurity_solvers.cthyb_krylov import *
from pytriqs.gf.local import *

# import parameters from cwd
from os import getcwd
from sys import path
path.insert(0,getcwd())
from params import *
del path[0]

def print_master(msg):
    if mpi.rank==0: print msg

pp = Parameters()
for k in p: pp[k] = p[k]

print_master("Welcome to spinless (spinless electrons on a correlated dimer) test.")

H = U*N("tot",0)*N("tot",1) -mu*(N("tot",0) + N("tot",1)) -t*(C_dag("tot",0)*C("tot",1) + C_dag("tot",1)*C("tot",0))
  
QN = []
if use_qn: QN.append(N("tot",0)+N("tot",1))

gf_struct = {"tot":[0,1]}

print_master("Constructing the solver...")

## Construct the solver
S = Solver(parameters=pp, H_local=H, quantum_numbers=QN, gf_struct=gf_struct)

print_master("Preparing the hybridization function...")

## Set hybridization function    
delta_w = GfImFreq(indices = [0,1], beta=beta)
delta_w <<= inverse(iOmega_n - numpy.matrix([[epsilon,-t],[-t,epsilon]])) + inverse(iOmega_n - numpy.matrix([[-epsilon,-t],[-t,-epsilon]]))
S.Delta_tau["tot"] <<= InverseFourier(delta_w)

print_master("Running the simulation...")

## Solve the problem
S.solve(parameters=pp)

## Save the results  
if mpi.rank==0:
    Results = HDFArchive(results_file_name(use_qn),'w')
    Results["tot"] = S.G_tau["tot"]