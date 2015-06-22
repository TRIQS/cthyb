#!/bin/env pytriqs

from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import *
from matplotlib.backends.backend_pdf import PdfPages
from itertools import product

def setup_fig():
    axes = plt.gca()
    axes.set_ylabel('$G(\\tau)$')
    axes.set_ylim(-1.0,.0)
    axes.legend(loc='lower center',prop={'size':10})

spin_names = ("up","dn")
num_orbitals = 2

pp = PdfPages('G.pdf')
ed_arch = HDFArchive('kanamori_offdiag.ed.h5','r')

for use_qn in (True,False):
    file_name = "kanamori_offdiag"
    if use_qn: file_name += ".qn"
    file_name += ".h5"

    try:
        arch = HDFArchive(file_name,'r')
        plt.clf()

        name = "cthyb (QN)" if use_qn else "cthyb"

        GF_up = rebinning_tau(arch['G_tau']['up'],200)
        GF_dn = rebinning_tau(arch['G_tau']['dn'],200)

        for o1, o2 in product(range(num_orbitals),range(num_orbitals)):
            oplot(GF_up[o1,o2], name=name+",$\uparrow%i%i$"%(o1,o2))
            oplot(GF_dn[o1,o2], name=name+",$\downarrow%i%i$"%(o1,o2))
            oplot(ed_arch['up'][o1,o2], name="ED,$\uparrow%i%i$"%(o1,o2))
            oplot(ed_arch['dn'][o1,o2], name="ED,$\downarrow%i%i$"%(o1,o2))

        setup_fig()
        pp.savefig(plt.gcf())

    except IOError: pass

pp.close()
