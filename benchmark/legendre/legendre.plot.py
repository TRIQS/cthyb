#!/bin/env python

import numpy as np
from h5 import *
from triqs.gf import *
from triqs.plot.mpl_interface import plt

arch = HDFArchive('legendre.h5','r')
ed_arch = HDFArchive('legendre.ed.h5','r')

plt.figure()
subp = [2, 1, 1]

plt.subplot(*subp); subp[-1] += 1

for spin in ("up","dn"):
    plt.plot(arch['G_l'][spin].data.flatten(),label="cthyb," + {'up':"$\uparrow\uparrow$",'dn':"$\downarrow\downarrow$"}[spin])
    plt.plot(ed_arch[spin].data.flatten(),label="ED," + {'up':"$\uparrow\uparrow$",'dn':"$\downarrow\downarrow$"}[spin])

axes = plt.gca()
axes.set_xlabel('$l$')
axes.set_ylabel('$G_l$')
axes.set_ylim((-2.0,1.0))
axes.legend(loc='lower center',prop={'size':10})

plt.subplot(*subp); subp[-1] += 1

for spin in ("up","dn"):
    plt.plot(np.abs(arch['G_l'][spin].data.flatten()),label="cthyb," + {'up':"$\uparrow\uparrow$",'dn':"$\downarrow\downarrow$"}[spin])
    plt.plot(np.abs(ed_arch[spin].data.flatten()),label="ED," + {'up':"$\uparrow\uparrow$",'dn':"$\downarrow\downarrow$"}[spin])

axes = plt.gca()
axes.set_xlabel('$l$')
axes.set_ylabel('$G_l$')
axes.legend(loc='best',prop={'size':10})
plt.semilogy([], [])

plt.tight_layout()
plt.savefig('G_l.pdf')
