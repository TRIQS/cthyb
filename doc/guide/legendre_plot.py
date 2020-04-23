from pytriqs.gf import *
from h5 import *
from pytriqs.plot.mpl_interface import oplot

A = HDFArchive("aim_solution.h5",'r')
oplot(A['G_l']['up'], '-o', x_window=(15,45) )
