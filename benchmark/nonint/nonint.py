#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.parameters.parameters import Parameters
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.operators import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import *
from matplotlib.backends.backend_pdf import PdfPages

def print_master(msg):
	if mpi.rank==0: print msg

print_master("Welcome to nonint (non-interacting many-band systems) test.")
print_master("This test is aimed to reveal excessive state truncation issues.")

beta = 40.0

for modes in range(1,10+1):
    V = [0.2]*modes
    e = [-0.2]*modes

    gf_struct = {str(n):[0] for n in range(0,len(V))}

    n_iw = 1025
    n_tau = 10001

    p = SolverCore.solve_parameters()
    p["max_time"] = -1
    p["random_name"] = ""
    p["random_seed"] = 123 * mpi.rank + 567
    p["verbosity"] = 2
    p["length_cycle"] = 50
    p["n_warmup_cycles"] = 50000
    p["n_cycles"] = 1200000

    # Local Hamiltonian
    H = Operator()

    # Quantum numbers (N_up and N_down)
    QN=[]
    for b in sorted(gf_struct.keys()): QN.append(n(b,0))

    print_master("Constructing the solver...")

    # Construct the solver
    S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

    print_master("Preparing the hybridization function...")

    # Set hybridization function
    for m, b in enumerate(sorted(gf_struct.keys())):
        delta_w = GfImFreq(indices = [0], beta=beta)
        delta_w <<= (V[m]**2) * inverse(iOmega_n - e[m])
        S.G0_iw[b][0,0] <<= inverse(iOmega_n - e[m] - delta_w)

    print_master("Running the simulation...")

    # Solve the problem
    S.solve(h_loc=H, params=p, quantum_numbers=QN, use_quantum_numbers=True)

    # Save and plot the results
    if mpi.rank==0:
        Results = HDFArchive('nonint_%i.h5'%len(V),'w')
        pdf = PdfPages("G_nonint_%i.pdf"%len(V))

        for m, b in enumerate(sorted(gf_struct.keys())):
            Results['G_'+b] = S.G_tau[b]
       
            g_theor = GfImTime(indices = [0], beta=beta, n_points=n_tau)
            e1 = e[m] - V[m]
            e2 = e[m] + V[m]
            g_theor_w = GfImFreq(indices = [0], beta=beta)
            g_theor_w <<= 0.5*inverse(iOmega_n - e1) + 0.5*inverse(iOmega_n - e2)
            g_theor[0,0] <<= InverseFourier(g_theor_w)
       
            plt.clf()
            oplot(S.G_tau[b][0,0], name="cthyb")
            oplot(g_theor[0,0], name="Theory")

        pdf.savefig(plt.gcf())

        pdf.close()
