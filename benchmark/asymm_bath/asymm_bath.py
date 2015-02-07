#!/bin/env pytriqs

import numpy as np
import pytriqs.utility.mpi as mpi
from pytriqs.gf.local import *
from pytriqs.operators import *
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb import *
import itertools

from pytriqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages

def print_master(msg):
    if mpi.rank==0: print msg

print_master("Welcome to asymm_bath test (1 band with a small asymmetric hybridization function).")
print_master("This test helps to detect sampling problems.")

# H_loc parameters
beta = 40.0
ed = -1.0
U = 2.0

epsilon = [0.0,0.1,0.2,0.3,0.4,0.5]
V = 0.2

# Parameters
n_iw = 1025
n_tau = 10001

p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 20000
p["n_cycles"] = 1000000
p["make_histograms"] = True

# Block structure of GF
gf_struct = {'up':[0], 'dn' : [0]}

# Hamiltonian
H = U*n("up",0)*n("dn",0)

# Quantum numbers
qn = []
qn.append(n("up",0))
qn.append(n("dn",0))
p["partition_method"] = "quantum_numbers"
p["quantum_numbers"] = qn

# Construct solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

def read_histo(f):
    histo = []
    for line in f:
        cols = filter(lambda s: s, line.split(' '))
        histo.append(float(cols[1]))
    for x in reversed(histo):
        if x != 0: break
        histo.remove(x)
    return histo

if mpi.rank==0:
    arch = HDFArchive('asymm_bath.h5','w')
    pp = PdfPages('G_asymm_bath.pdf')

# Set hybridization function
for e in epsilon:
    delta_w = GfImFreq(indices = [0], beta=beta)
    delta_w << (V**2) * inverse(iOmega_n - e)

    S.G0_iw["up"] << inverse(iOmega_n - ed - delta_w)
    S.G0_iw["dn"] << inverse(iOmega_n - ed - delta_w)

    S.solve(h_loc=H, **p)

    if mpi.rank==0:
        arch['epsilon_' + str(e)] = {"up":S.G_tau["up"], "dn":S.G_tau["dn"]}

        plt.clf()
        oplot(rebinning_tau(S.G_tau['up'],300), name="$\uparrow\uparrow$")
        oplot(rebinning_tau(S.G_tau['dn'],300),name="$\downarrow\downarrow$")

        a = plt.gca()
        a.set_ylabel('$G(\\tau)$')
        a.set_xlim((0,beta))
        a.set_ylim((-1,0))
        a.legend(loc='lower right',prop={'size':10})

        a.set_title("$U=%.1f$, $\epsilon_d=%.1f$, $V=%.1f$, $\epsilon_k=%.1f$" % (U,ed,V,e))

        histo = read_histo(open("histo_opcount_total.dat",'r'))

        histo_a = plt.axes([.35, .15, .3, .4], axisbg='y')
        histo_a.bar(range(len(histo)), histo)
        histo_a.set_title('histo_opcount_total')

        pp.savefig(plt.gcf())

if mpi.rank==0: pp.close()
