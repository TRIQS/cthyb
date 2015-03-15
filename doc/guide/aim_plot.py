from pytriqs.gf.local import *
from pytriqs.archive import *
from pytriqs.plot.mpl_interface import oplot

A = HDFArchive('aim_solution.h5','r')
oplot(A['G_iw']['up'], '-o', x_window = (0,10))
