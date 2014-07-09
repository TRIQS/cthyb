#!/bin/env pytriqs

from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import *
from matplotlib.backends.backend_pdf import PdfPages

def setup_fig():
    axes = plt.gca()
    axes.set_ylabel('$G(\\tau)$')
    axes.set_ylim(-1.0,.0)
    axes.legend(loc='lower center',prop={'size':10})

matrix={   
    'file':"kanamori.matrix.h5",
        'name':"matrix",
        'path':["%(spin)s-%(orbital)i"],
}
ed={    
    'file':"kanamori.ed.h5",
    'name':"ED",
    'path':["%(spin)s-%(orbital)i"],
}
cthyb_qn={
    'file':"kanamori.qn.h5",
    'name':"cthyb (QN)",
    'path':["%(spin)s-%(orbital)i"]
}
cthyb={
    'file':"kanamori.h5",
    'name':"cthyb",
    'path':["%(spin)s-%(orbital)i"]
}

objects_to_plot = [[matrix,cthyb],[matrix,cthyb_qn],[ed,cthyb],[ed,cthyb_qn],[matrix,ed]]

num_o = 2

pp = PdfPages('G.pdf')
for plot_objs in objects_to_plot:
    try:
        plt.clf()
        for obj in plot_objs:
            arch = HDFArchive(obj['file'],'r')
            
            for o in range(0,num_o):
                target_up = arch
                for p in obj['path']: target_up = target_up[p % {'spin':'up','orbital':o}]    
                oplot(target_up, name=obj['name'] + (",$\uparrow%i$"%o))
                target_down = arch
                for p in obj['path']: target_down = target_down[p % {'spin':'dn','orbital':o}]    
                oplot(target_down, name=obj['name'] + (",$\downarrow%i$"%o))
                        
        setup_fig()
        pp.savefig(plt.gcf())

    except IOError: pass

pp.close()