from pytriqs.gf import *
from h5 import *
from pytriqs.plot.mpl_interface import oplot, plt

with HDFArchive("dmft_solution.h5",'r') as ar:
    for i in range(5):
        oplot(ar['G_iw-%i'%i]['up'], '-o', mode = 'I',
              label = 'Iteration = %s'%i,
              x_window = (0,2))

plt.legend(loc=4)
