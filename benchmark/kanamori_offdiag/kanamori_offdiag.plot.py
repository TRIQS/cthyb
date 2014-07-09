#!/bin/env pytriqs

from itertools import product
from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import *
from matplotlib.backends.backend_pdf import PdfPages

def setup_fig():
    axes = plt.gca()
    axes.set_ylabel('$G(\\tau)$')
    axes.set_ylim(-1.0,1.0)
    axes.legend(loc='lower center',prop={'size':10},ncol=2)

matrix={   
    'file':"kanamori_offdiag.matrix.h5",
    'name':"matrix",
    'path':["%(spin)s"],
}
ed={    
    'file':"kanamori_offdiag.ed.h5",
    'name':"ED",
    'path':["%(spin)s"],
}
cthyb_qn={
    'file':"kanamori_offdiag.qn.h5",
    'name':"cthyb (QN)",
    'path':["%(spin)s"]
}
cthyb={
    'file':"kanamori_offdiag.h5",
    'name':"cthyb",
    'path':["%(spin)s"]
}

objects_to_plot = [[matrix,cthyb],[matrix,cthyb_qn],[cthyb,ed],[cthyb_qn,ed],[matrix,ed]]

num_o = 2

pp = PdfPages('G.pdf')
for plot_objs in objects_to_plot:
    try:
        plt.clf()
        for obj in plot_objs:
            arch = HDFArchive(obj['file'],'r')

            target_up = arch
            for p in obj['path']: target_up = target_up[p % {'spin':'up'}]    
            for o1,o2 in product(range(0,num_o),range(0,num_o)):
                oplot(target_up[o1,o2], name=obj['name'] + (",$\uparrow%i%i$"%(o1,o2)))
            
            target_down = arch
            for p in obj['path']: target_down = target_down[p % {'spin':'dn'}]
            for o1,o2 in product(range(0,num_o),range(0,num_o)):
                oplot(target_down[o1,o2], name=obj['name'] + (",$\downarrow%i%i$"%(o1,o2)))
            
        setup_fig()
        pp.savefig(plt.gcf())

    except IOError: pass
    
pp.close()
