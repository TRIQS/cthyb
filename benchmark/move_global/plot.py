#!/bin/env pytriqs

from pytriqs.archive import HDFArchive
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from sys import argv

results_file = argv[1]

float_cmp = lambda x,y: cmp(float(x), float(y))

arch = HDFArchive(results_file,'r')
beta = arch['beta']
move_global_prob = arch['move_global_prob']

all_mg = arch.keys()
all_mg.remove('beta')
all_mg.remove('move_global_prob')

pp = PdfPages('occupation_mg_beta%.0f_prob%.2f.pdf' % (beta,move_global_prob))

dat_N1_up = {mg: arch[mg]['N1_up'] for mg in all_mg}
dat_N1_dn = {mg: arch[mg]['N1_dn'] for mg in all_mg}
dat_N2_up = {mg: arch[mg]['N2_up'] for mg in all_mg}
dat_N2_dn = {mg: arch[mg]['N2_dn'] for mg in all_mg}

x = sorted(dat_N1_up.keys())
index = np.arange(len(x))
bar_width = 0.5

y_N1_dn = [dat_N1_dn[str(key)] for key in x]
y_N1_up = [dat_N1_up[str(key)] for key in x]
y_N2_dn = [dat_N2_dn[str(key)] for key in x]
y_N2_up = [dat_N2_up[str(key)] for key in x]

plt.title("Test move_global $\\beta$ = %.0f, move_global_prob = %.3f" % (beta,move_global_prob), fontsize=12)
plt.xlim(-0.35,4.35)
plt.ylim(0,1.2)
plt.xticks(index, x)
plt.axhline(0.5, color='grey', linewidth = 0.25, linestyle = '--')
plt.axhline(1.0, color='grey')
plt.bar(index-0.13, y_N1_dn, align='center', width=0.15, color = 'royalblue', label = '$\\downarrow$')
plt.bar(index+0.13, y_N2_dn, align='center', width=0.15, color = 'royalblue')
plt.bar(index-0.13, y_N1_up, bottom=y_N1_dn, align='center', width=0.15, color = 'tomato', label = '$\\uparrow$' )
plt.bar(index+0.13, y_N2_up, bottom=y_N2_dn, align='center', width=0.15, color = 'tomato')
plt.legend(loc = 'upper center', fontsize = 9)
pp.savefig()

pp.close()

