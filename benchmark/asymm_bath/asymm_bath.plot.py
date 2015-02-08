#!/bin/env pytriqs

from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
import itertools

from pytriqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages

arch = HDFArchive('asymm_bath.h5','r')
arch_ed = HDFArchive('asymm_bath.ed.h5','r')

# Plot G(\tau)
pp = PdfPages('G_asymm_bath.pdf')

for e_group_name in arch:
    e_group = arch[e_group_name]
    e_group_ed = arch_ed[e_group_name]

    beta = e_group['beta']
    U = e_group['U']
    ed = e_group['ed']
    V = e_group['V']
    e = e_group['e']

    plt.clf()
    oplot(rebinning_tau(e_group['G_tau']['up'],300),name="CTHYB,$\uparrow\uparrow$")
    oplot(rebinning_tau(e_group['G_tau']['dn'],300),name="CTHYB,$\downarrow\downarrow$")

    oplot(rebinning_tau(e_group_ed['G_tau']['up'],300),name="ED,$\uparrow\uparrow$")
    oplot(rebinning_tau(e_group_ed['G_tau']['dn'],300),name="ED,$\downarrow\downarrow$")

    a = plt.gca()
    a.set_ylabel('$G(\\tau)$')
    a.set_xlim((0,beta))
    a.set_ylim((-1,0))
    a.legend(loc='lower right',prop={'size':10})

    a.set_title("$U=%.1f$, $\epsilon_d=%.1f$, $V=%.1f$, $\epsilon_k=%.1f$" % (U,ed,V,e))

    histo_a = plt.axes([.35, .15, .3, .4], axisbg='y')
    opcount_data = []
    for opcount, w, int_w in reversed(e_group['histo_opcount_total']):
        if w==0: continue
        opcount_data.insert(0,w)
    histo_a.bar(range(len(opcount_data)), opcount_data)
    histo_a.set_title('histo_opcount_total')

    pp.savefig(plt.gcf())

pp.close()

# Plot length histograms
pp = PdfPages('length_histos_asymm_bath.pdf')

for e_group_name in arch:
    e_group = arch[e_group_name]

    beta = e_group['beta']
    U = e_group['U']
    ed = e_group['ed']
    V = e_group['V']
    e = e_group['e']

    def plot_histos(title,histos):
        a = plt.gca()
        a.set_xlim((0,beta))
        a.set_xlabel('$\\Delta\\tau$')
        a.set_ylabel(title)

        for name, histo in histos.items():
            dtau, w, int_w = zip(*histo)
            plt.plot(dtau,w,label=name,linewidth=0.7)

        a.legend(loc='upper center',prop={'size':10})

    plt.clf()

    plt.suptitle("$U=%.1f$, $\epsilon_d=%.1f$, $V=%.1f$, $\epsilon_k=%.1f$" % (U,ed,V,e))

    # Move insert
    plt.subplot(3,1,1)
    plot_histos("Insertion",{"Proposed" : e_group['histo_insert_length_proposed'],
                             "Accepted" : e_group['histo_insert_length_accepted']})
    # Move remove
    plt.subplot(3,1,2)
    plot_histos("Removal",{"Proposed" : e_group['histo_remove_length_proposed'],
                           "Accepted" : e_group['histo_remove_length_accepted']})
    # Move shift
    plt.subplot(3,1,3)
    plot_histos("Shift",{"Proposed" : e_group['histo_shift_length_proposed'],
                         "Accepted" : e_group['histo_shift_length_accepted']})

    pp.savefig(plt.gcf())

pp.close()