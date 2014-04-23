#!/bin/env pytriqs

from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages

# import parameters from cwd
from os import getcwd
from sys import path
path.insert(0,getcwd())
import params
del path[0]

objects_to_plot=[['matrix',(False,False)],['matrix',(False,True)],['matrix',(True,False)],['matrix',(True,True)],
                 ['ed',(False,False)],['ed',(False,True)],['ed',(True,False)],['ed',(True,True)],
                 ['ed','matrix']]

def setup_fig():
    axes = plt.gca()
    axes.set_ylabel('$G(\\tau)$')
    axes.legend(loc='lower center',prop={'size':10})

pp = PdfPages('G.pdf')
for plot_objs in objects_to_plot:
    try:
        plt.clf()
        for obj in plot_objs:

            if type(obj) is tuple:
                filename = params.results_file_name(*obj)
                params.use_blocks = obj[0]
                
                name_parts = []
                if obj[0]: name_parts.append('Block')
                if obj[1]: name_parts.append('QN')
                name = 'cthyb' + (' (' + ', '.join(name_parts) + ')' if len(name_parts) else '')
                
            else:
                filename = 'anderson.' + obj + '.h5'
                name = {'ed':'ED', 'matrix':'Matrix'}[obj]
                params.use_blocks = True
            
            arch = HDFArchive(filename,'r')
            
            for spin in params.spin_names:
                bn, i = params.mkind(spin)
                oplot(arch[bn][i,i], name=name + "," + {'up':"$\uparrow\uparrow$",'dn':"$\downarrow\downarrow$"}[spin])
                            
        setup_fig()
        pp.savefig(plt.gcf())

    except IOError: pass

pp.close()
