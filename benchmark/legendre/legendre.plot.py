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

objects_to_plot=[('cthyb','matrix'),('cthyb','ED'),('matrix','ED')]

def setup_fig():
    axes = plt.gca()
    axes.set_xlabel('$l$')
    axes.set_ylabel('$G_l$')
    axes.set_ylim((-2.0,1.0))
    axes.legend(loc='lower center',prop={'size':10})

pp = PdfPages('G.pdf')
for plot_objs in objects_to_plot:
    try:
        plt.clf()
        for obj in plot_objs:

            filename = {'cthyb' : params.results_file_name,
                        'matrix' : params.matrix_results_file_name,
                        'ED' : params.ref_file_name }[obj]

            arch = HDFArchive(filename,'r')

            for spin in params.spin_names:
                bn, i = params.mkind(spin)
                plt.plot(arch[bn].data.flatten(),label=obj + "," + {'up':"$\uparrow\uparrow$",'dn':"$\downarrow\downarrow$"}[spin])

        setup_fig()
        pp.savefig(plt.gcf())

    except IOError: pass

pp.close()
