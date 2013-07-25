from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import *

plot_num = 3

# Anderson test
if plot_num == 1:
  A = HDFArchive("anderson.krylov.h5",'r')
  B = HDFArchive("anderson.hyb.h5",'r')
  oplot(A['G_up'])
  oplot(A['G_down'])
  oplot(B['G_tau'])

# Spinless test
if plot_num == 2:
  A = HDFArchive("spinless.krylov.h5",'r')
  B = HDFArchive("spinless.hyb.h5",'r')
  oplot(A['G_tau'])
  oplot(B['G_tau'])

# Kanamori test
if plot_num == 3:
  A = HDFArchive("kanamori.krylov.h5",'r')
  B = HDFArchive("kanamori.hyb.h5",'r')
  oplot(A['G_up-0'])
  oplot(A['G_down-0'])
  oplot(A['G_up-1'])
  oplot(A['G_down-1'])
  oplot(B['G_tau'])
