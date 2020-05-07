#!/bin/env python

from h5 import *
from triqs.gf import *
from triqs.gf.gf_fnt import rebinning_tau
from triqs.plot.mpl_interface import *
from matplotlib.backends.backend_pdf import PdfPages
from itertools import product

def setup_fig():
    axes = plt.gca()
    axes.set_ylabel('$G(\\tau)$')
    axes.legend(loc='lower center',prop={'size':8}, ncol=1)

spin_names = ("up","dn")
num_orbitals = 2

pp = PdfPages('G.pdf')
ed_arch = HDFArchive('kanamori_offdiag.ed.h5','r')
pyed_arch = HDFArchive('kanamori_offdiag.pyed.h5','r')
pomerol_arch = HDFArchive('kanamori_offdiag.pomerol.h5','r')

for use_qn in (True,False):
    file_name = "kanamori_offdiag"
    if use_qn: file_name += ".qn"
    file_name += ".h5"

    try:
        arch = HDFArchive(file_name,'r')

        name = "cthyb (QN)" if use_qn else "cthyb"

        GF_up = rebinning_tau(arch['G_tau']['up'],200)
        GF_dn = rebinning_tau(arch['G_tau']['dn'],200)

        ed_opt = dict(lw=2.0, alpha=1.0)
        cthyb_opt = dict(lw=1.0, alpha=1.0)
        
        for o1, o2 in product(range(num_orbitals), repeat=2):
            plt.clf()
            plt.title('using_qn = ' + str(use_qn))
            oplot(ed_arch['up'][o1,o2], name="ED,$\uparrow%i%i$"%(o1,o2), **ed_opt)
            oplot(pyed_arch['up'][o1,o2], name="PYED,$\uparrow%i%i$"%(o1,o2), **ed_opt)
            oplot(pomerol_arch['up']['up'][o1,o2], 'o', name="Pomerol,$\uparrow%i%i$"%(o1,o2), **ed_opt)
            oplotr(GF_up[o1,o2], name=name+",$\uparrow%i%i$"%(o1,o2), **cthyb_opt)
            oploti(GF_up[o1,o2], name=name+",$\uparrow%i%i$"%(o1,o2), **cthyb_opt)
            setup_fig()
            pp.savefig(plt.gcf())

            plt.clf()
            plt.title('using_qn = ' + str(use_qn))
            oplot(ed_arch['dn'][o1,o2], name="ED,$\downarrow%i%i$"%(o1,o2), **ed_opt)
            oplot(pyed_arch['dn'][o1,o2], name="PYED,$\downarrow%i%i$"%(o1,o2), **ed_opt)
            oplot(pomerol_arch['dn']['do'][o1,o2], 'o', name="Pomerol,$\uparrow%i%i$"%(o1,o2), **ed_opt)
            oplotr(GF_dn[o1,o2], name=name+",$\downarrow%i%i$"%(o1,o2), **cthyb_opt)
            oploti(GF_dn[o1,o2], name=name+",$\downarrow%i%i$"%(o1,o2), **cthyb_opt)
            setup_fig()
            pp.savefig(plt.gcf())
            
    except IOError: pass

pp.close()
