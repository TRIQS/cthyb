#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi
from pytriqs.archive import HDFArchive
from pytriqs.parameters.parameters import Parameters
from pytriqs.applications.impurity_solvers.cthyb_krylov import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import *
from matplotlib.backends.backend_pdf import PdfPages

def print_master(msg):
	if mpi.rank==0: print msg

beta = 40.0

for modes in range(1,10+1):
	V = [0.2]*modes
	e = [-0.2]*modes

	gf_struct = {str(n):[0] for n in range(0,len(V))}

	p = {}
	p["beta"] = beta
	p["max_time"] = -1
	p["random_name"] = ""
	p["random_seed"] = 123 * mpi.rank + 567
	p["verbosity"] = 3
	p["length_cycle"] = 50
	p["n_warmup_cycles"] = 10000
	p["n_cycles"] = 300000
	p["n_tau_delta"] = 1000
	p["n_tau_g"] = 400
	p["krylov_bs_use_cutoff"] = True
	p["krylov_bs_prob_cutoff"] = -1.0
	p["krylov_gs_energy_convergence"] = 1e-8
	p["krylov_small_matrix_size"] = 100

	pp = Parameters()
	for k in p: pp[k] = p[k]

	# Local Hamiltonian
	H = Operator()
	for n, b in enumerate(sorted(gf_struct.keys())):
    		H += e[n]*N(b,0)

	# Quantum numbers (N_up and N_down)
	QN=[]
	for b in sorted(gf_struct.keys()): QN.append(N(b,0))

	print_master("Constructing the solver...")

	# Construct the solver
	S = Solver(parameters=pp, H_local=H, quantum_numbers=QN, gf_struct=gf_struct)

	print_master("Preparing the hybridization function...")

	# Set hybridization function
    
	for n, b in enumerate(sorted(gf_struct.keys())):
    		delta_w = GfImFreq(indices = [0], beta=beta)
    		delta_w <<= (V[n]**2) * inverse(iOmega_n - e[n])
    		S.Delta_tau[b][0,0] <<= InverseFourier(delta_w)

	print_master("Running the simulation...")

	# Solve the problem
	S.solve(parameters=pp)

	# Save and plot the results  
	if mpi.rank==0:
    		Results = HDFArchive('nonint_%i.h5'%len(V),'w')
    		pdf = PdfPages("G_nonint_%i.pdf"%len(V))
    
    		for n, b in enumerate(sorted(gf_struct.keys())):
    			Results['G_'+b] = S.G_tau[b]

    			g_theor = GfImTime(indices = [0], beta=beta, n_points=p["n_tau_g"])
    			e1 = e[n] - V[n]
    			e2 = e[n] + V[n]
    			g_theor_w = GfImFreq(indices = [0], beta=beta)
    			g_theor_w <<= 0.5*inverse(iOmega_n - e1) + 0.5*inverse(iOmega_n - e2)
    			g_theor[0,0] <<= InverseFourier(g_theor_w)

    			plt.clf()
			oplot(S.G_tau[b][0,0], name="Krylov")
    			oplot(g_theor[0,0], name="Theory")
    
    			pdf.savefig(plt.gcf())
    
    		pdf.close()
