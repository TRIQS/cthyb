from pytriqs.gf.local import *
from pytriqs.archive import *
from pytriqs.plot.mpl_interface import oplot

A = HDFArchive("aim_solution.h5")
oplot(A['G_l']['up'], '-o', x_window=(15,45) )
