#!/bin/env python

from h5 import HDFArchive
from triqs.gf import *
from triqs.gf.gf_fnt import rebinning_tau
from triqs.plot.mpl_interface import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from sys import argv
from itertools import product

plot_n_tau = 400

results_file = argv[1]

float_cmp = lambda x,y: cmp(float(x), float(y))

arch = HDFArchive(results_file,'r')
beta = arch['beta']
move_global_prob = arch['move_global_prob']

all_mg = list(arch.keys())
all_mg.remove('beta')
all_mg.remove('move_global_prob')

# Plot occupations
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

# Plot Green's functions
pp = PdfPages('G_tau_mg_beta%.0f_prob%.2f.pdf' % (beta,move_global_prob))

mg_pairs = [(a,b) for a,b in product(all_mg,all_mg) if a<b]
spin_labels = {'up':'$\uparrow\uparrow$', 'dn':'$\downarrow\downarrow$'}
for mg1, mg2 in mg_pairs:

    plt.clf()
    for s, i in product(('up','dn'),(0,1)):
        G_tau_1 = rebinning_tau(arch[mg1]['G_tau'][s], plot_n_tau)
        G_tau_2 = rebinning_tau(arch[mg2]['G_tau'][s], plot_n_tau)

        oplot(G_tau_1[i,i], alpha=0.25, mode='R', label="%s,%s,%i%i" % (mg1,s,i,i))
        oplot(G_tau_2[i,i], alpha=0.25, mode='R', label="%s,%s,%i%i" % (mg2,s,i,i))

    plt.title("Test move_global $\\beta$ = %.0f, move_global_prob = %.3f" % (beta,move_global_prob), fontsize=12)
    plt.ylim(-0.6,0.05)
    plt.ylabel('$G(\\tau)$')
    plt.legend(loc = 'lower center', fontsize = 9)
    pp.savefig()

pp.close()

