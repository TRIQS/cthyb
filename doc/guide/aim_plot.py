from pytriqs.gf import *
from pytriqs.archive import *
from pytriqs.plot.mpl_interface import oplot

with HDFArchive('aim_solution.h5','r') as ar:
    oplot(ar['G_iw']['up'], '-o', x_window = (0,10))
