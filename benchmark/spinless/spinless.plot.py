#!/bin/env pytriqs

from itertools import product
from pytriqs.archive import *
from pytriqs.gf import *
from pytriqs.gf.gf_fnt import rebinning_tau
from pytriqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages

spin_names = ("up","dn")

def setup_fig():
    axes = plt.gca()
    axes.set_ylabel('$G(\\tau)$')
    axes.legend(loc='best',prop={'size':8})

pp = PdfPages('G.pdf')
ed_arch = HDFArchive('spinless.ed.h5','r')

for use_qn in (False,True):
    file_name = "spinless"
    if use_qn: file_name += ".qn"
    file_name += ".h5"

    try:
        arch = HDFArchive(file_name,'r')
        plt.clf()

        name = 'cthyb' + (' (QN)' if use_qn else '')
        for i1,i2 in product(("A","B"),("A","B")):
            GF = rebinning_tau(arch['tot'][i1,i2],500)
            oplot(arch["tot"][i1,i2], name=name + ",%s%s" % (i1,i2))
            oplot(ed_arch["tot"][i1,i2], name="ED, %s%s" % (i1,i2))

        setup_fig()
        pp.savefig(plt.gcf())

    except IOError: pass

pp.close()
