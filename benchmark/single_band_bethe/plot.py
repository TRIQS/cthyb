
import numpy as np
from pytriqs.gf import *
from pytriqs.archive import HDFArchive
from pytriqs.plot.mpl_interface import oplot, oplotr, oploti, plt

def plot_full_iter(iter_no, filename):

    i = iter_no
    
    with HDFArchive(filename, 'r') as a:

        U = a['U']

        G0_iw = a['G0_iw-%i' % i]
        Delta_tau = a['Delta_tau-%i' % i]

        G_iw = a['G_iw-%i' % i]
        G_iw_raw = a['G_iw_raw-%i' % i]

        Sigma_iw = a['Sigma_iw-%i' % i]
        Sigma_iw_raw = a['Sigma_iw_raw-%i' % i]


    Delta_0_iw = G_iw.copy()
    for name, d in Delta_0_iw:
        d << inverse(G0_iw[name]) - iOmega_n - 0.5*U

    Delta_1_iw = G_iw.copy()
    for name, d in Delta_1_iw:
        d << inverse(G_iw[name]) - iOmega_n - 0.5*U + Sigma_iw[name]

    #Delta_iw_ref = G_iw.copy()
    #for name, d in Delta_iw_ref:
    #    d << Fourier(Delta_tau[name])

    plt.figure(figsize=(10, 10))
    subp = [4, 4, 1]

    plt.subplot(*subp); subp[-1] += 1
    plt.title('Delta_tau')
    oplot(Delta_tau) 
    ax = plt.gca()
    ax.legend_ = None

    plt.subplot(*subp); subp[-1] += 1
    plt.title('G0_iw')
    oplot(G0_iw)
    ax = plt.gca()
    ax.legend_ = None

    plt.subplot(*subp); subp[-1] += 1
    plt.title('G_iw_raw')
    oplot(G_iw_raw)
    ax = plt.gca()
    ax.legend_ = None

    plt.subplot(*subp); subp[-1] += 1
    plt.title('G_iw')
    oplot(G_iw)
    ax = plt.gca()
    ax.legend_ = None

    plt.subplot(*subp); subp[-1] += 1
    plt.title('Sigma_iw_raw')
    oplot(Sigma_iw_raw)
    ax = plt.gca()
    ax.legend_ = None

    plt.subplot(*subp); subp[-1] += 1
    plt.title('Sigma_iw_raw')
    oplot(Sigma_iw_raw)
    ax = plt.gca()
    ax.legend_ = None
    plt.ylim([-6, 6])
    plt.xlim([-15, 15])
    
    plt.subplot(*subp); subp[-1] += 1
    plt.title('Sigma_iw')
    oplot(Sigma_iw)
    ax = plt.gca()
    ax.legend_ = None
    
    plt.subplot(*subp); subp[-1] += 1
    plt.title('Delta_0_iw, Delta_1_iw')
    oplot(Delta_0_iw)
    oplot(Delta_1_iw)
    ax = plt.gca()
    ax.legend_ = None

    plt.subplot(*subp); subp[-1] += 1
    plt.title('Delta_0_iw - Delta_1_iw')
    oplot(Delta_0_iw - Delta_1_iw)
    ax = plt.gca()
    ax.legend_ = None

    plt.subplot(*subp); subp[-1] += 1
    plt.title('Delta_0_iw + G_iw')
    oplot(Delta_0_iw + G_iw)
    ax = plt.gca()
    ax.legend_ = None

    plt.tight_layout()

def plot_all_iter(iter_min=0, iter_max=0, filename=None):

    plt.figure(figsize=(10, 10))

    with HDFArchive(filename, 'r') as a:
        U = a['U']
        for i in xrange(iter_min, iter_max):

            Delta_tau = a['Delta_tau-%i' % i]
            oplotr(Delta_tau['up'], label=i, alpha=0.25)
            oplotr(Delta_tau['down'], label=i, alpha=0.25)

    plt.tight_layout()

    plt.figure(figsize=(10, 10))

    with HDFArchive(filename, 'r') as a:
        U = a['U']
        for i in xrange(iter_min, iter_max):
            G_iw = a['G_iw-%i' % i]
            oplotr(G_iw['up'], label=i, alpha=0.5)
            oploti(G_iw['up'], label=i, alpha=0.5)

    plt.tight_layout()

    plt.figure(figsize=(10, 10))

    with HDFArchive(filename, 'r') as a:
        U = a['U']
        for i in xrange(iter_min, iter_max):
            Sigma_iw = a['Sigma_iw-%i' % i]
            oplotr(Sigma_iw['up'], label=i, alpha=0.5)
            oploti(Sigma_iw['up'], label=i, alpha=0.5)

    plt.tight_layout()
    
#filename = 'data-U5.00_tail_fit_False.h5'
filename = 'data-U5.00_tail_fit_True.h5'
plot_full_iter(iter_no=8, filename=filename)
plot_all_iter(iter_min=3, iter_max=6, filename=filename)

plt.show()
