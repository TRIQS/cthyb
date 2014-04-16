#!/bin/env pytriqs

from itertools import product
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

objects_to_plot=[['matrix',(False,)],['matrix',(True,)],
                 ['ed',(False,)],['ed',(True,)],
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
                name = 'cthyb'
                if obj[0]: name += '(QN)'
            else:
                filename = 'spinless.' + obj + '.h5'
                name = {'ed':'ED', 'matrix':'Matrix'}[obj]
            
            arch = HDFArchive(filename,'r')
            
            for i1,i2 in product((0,1),(0,1)):
                oplot(arch["tot"][i1,i2], name=name + ",%i%i" % (i1,i2))
                            
        setup_fig()
        pp.savefig(plt.gcf())

    except IOError: pass

pp.close()
