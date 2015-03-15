from pytriqs.gf.local import *
from pytriqs.archive import *
from pytriqs.plot.mpl_interface import *

A = HDFArchive("dmft_solution.h5",'r')

for i in range(5): 
    oplot(A['G_iw-%s'%i]['up'].imag,'-o', name='Iteration = %s'%i, x_window = (0,2))
plt.legend(loc=4)
