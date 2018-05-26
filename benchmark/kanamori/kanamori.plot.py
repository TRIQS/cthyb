#!/bin/env pytriqs

from pytriqs.archive import *
from pytriqs.gf import *
from pytriqs.gf.gf_fnt import rebinning_tau
from pytriqs.plot.mpl_interface import *
from matplotlib.backends.backend_pdf import PdfPages

def setup_fig():
    axes = plt.gca()
    axes.set_ylabel('$G(\\tau)$')
    axes.set_ylim(-1.0,.1)
    axes.legend(loc='best',prop={'size':8},ncol=2)

spin_names = ("up","dn")
num_orbitals = 2

pp = PdfPages('G.pdf')
ed_arch = HDFArchive('kanamori.ed.h5','r')

for use_qn in (True,False):
    file_name = "kanamori"
    if use_qn: file_name += ".qn"
    file_name += ".h5"

    try:
        arch = HDFArchive(file_name,'r')
        plt.clf()

        name = "cthyb (QN)" if use_qn else "cthyb"
        for o in range(num_orbitals):
            oplot(rebinning_tau(arch['G_tau']['up_%i'%o],200), name=name+",$\uparrow%i$"%o)
            oplot(rebinning_tau(arch['G_tau']['dn_%i'%o],200), name=name+",$\downarrow%i$"%o)
            oplot(ed_arch['up-%i'%o], name="ED,$\uparrow%i$"%o)
            oplot(ed_arch['dn-%i'%o], name="ED,$\downarrow%i$"%o)

        setup_fig()
        pp.savefig(plt.gcf())

    except IOError: pass

pp.close()
