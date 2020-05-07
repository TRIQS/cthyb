#!/bin/env python

from h5 import *
from triqs.gf import *
from triqs.gf.gf_fnt import rebinning_tau
from triqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages

def setup_fig():
    axes = plt.gca()
    axes.set_ylabel('$G(\\tau)$')
    axes.legend(loc='lower center',prop={'size':10})

spin_names = ("up","dn")

pp = PdfPages('G.pdf')
ed_arch = HDFArchive('anderson.ed.h5','r')
pyed_arch = HDFArchive('anderson.pyed.h5','r')

for use_blocks, use_qn in ((False,False),(True,False),(False,True),(True,True)):
    file_name = "anderson"
    if use_blocks: file_name += ".block"
    if use_qn: file_name += ".qn"
    file_name += ".h5"

    mkind = lambda spin: (spin,0) if use_blocks else ("tot",spin)

    try:
        arch = HDFArchive(file_name,'r')
        plt.clf()

        name_parts = []
        if use_blocks: name_parts.append('Block')
        if use_qn: name_parts.append('QN')
        name = 'cthyb' + (' (' + ', '.join(name_parts) + ')' if len(name_parts) else '')

        for spin in spin_names:
            bn, i = mkind(spin)
            GF = rebinning_tau(arch['G_tau'][bn],500)
            if use_blocks:
                oplot(GF, name=name + "," + {'up':"$\\uparrow\\uparrow$",'dn':"$\downarrow\downarrow$"}[spin])
            else:
                i = spin_names.index(i)
                oplot(GF[i,i], name=name + "," + {'up':"$\\uparrow\\uparrow$",'dn':"$\downarrow\downarrow$"}[spin])
            oplot(ed_arch[spin], name="ED," + {'up':"$\\uparrow\\uparrow$",'dn':"$\downarrow\downarrow$"}[spin])
            oplot(pyed_arch[spin], name="PYED," + {'up':"$\\uparrow\\uparrow$",'dn':"$\downarrow\downarrow$"}[spin])

        setup_fig()
        pp.savefig(plt.gcf())

    except IOError: pass

pp.close()
