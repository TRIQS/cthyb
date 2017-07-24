#!/bin/env pytriqs

from pytriqs.archive import *
from pytriqs.gf import *
from pytriqs.gf.gf_fnt import rebinning_tau
from pytriqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages

N_max = 10
arch = HDFArchive('nonint.h5','r')

for i in arch:
    subarch = arch[i]
    pp = PdfPages("G_nonint_%s.pdf"%i)

    G_tau = subarch['G_tau']
    V = subarch['V']
    e = subarch['e']
    beta = G_tau.mesh.beta
    n_tau = len(G_tau.mesh)

    for m, b in enumerate(G_tau.indices):
        g_theor = GfImTime(indices = [0], beta=beta, n_points=n_tau)
        e1 = e[m] - V[m]
        e2 = e[m] + V[m]
        g_theor_w = GfImFreq(indices = [0], beta=beta)
        g_theor_w << 0.5*inverse(iOmega_n - e1) + 0.5*inverse(iOmega_n - e2)
        g_theor[0,0] << InverseFourier(g_theor_w)

        plt.clf()
        oplot(rebinning_tau(G_tau[b][0,0],200), name="cthyb")
        oplot(g_theor[0,0], name="Theory")

        pp.savefig(plt.gcf())

    pp.close()
