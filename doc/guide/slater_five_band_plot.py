from pytriqs.gf import *
from pytriqs.gf.gf_fnt import rebinning_tau
from pytriqs.archive import *
from pytriqs.plot.mpl_interface import oplot

with HDFArchive('slater_five_band.h5','r') as ar:
    # Calculate orbital- and spin-averaged Green's function
    G_tau = ar['G_tau-1']
    g_tau_ave = G_tau['up_0'].copy()
    g_tau_ave.zero()
    for name, g in G_tau: g_tau_ave += g
    g_tau_ave = g_tau_ave/10.
    g_tau_rebin = rebinning_tau(g_tau_ave,1000)
    g_tau_rebin.name = r'$G_{\rm ave}$'
    oplot(g_tau_rebin,linewidth=2,label='')
