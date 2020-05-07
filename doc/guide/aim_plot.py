from triqs.gf import *
from h5 import *
from triqs.plot.mpl_interface import oplot

with HDFArchive('aim_solution.h5','r') as ar:
    oplot(ar['G_iw']['up'], '-o', x_window = (0,10))
