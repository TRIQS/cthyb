@preamble@
from itertools import product
from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import *
from matplotlib.backends.backend_pdf import PdfPages

def setup_fig():
	axes = plt.gca()
	axes.set_ylabel('$G(\\tau)$')
	axes.legend(loc='lower center',prop={'size':10},ncol=2)

hyb={	
	'file':"kanamori_offdiag.hyb.h5",
     	'name':"HYB",
     	'path':["G_tau","%(spin)s"],
}
ed={	
	'file':"kanamori_offdiag.ed.h5",
	'name':"ED",
	'path':["G_%(spin)s"],
}
krylov_qn={
	'file':"kanamori_offdiag.krylov_qn.h5",
	'name':"Krylov (QN)",
	'path':["G_%(spin)s"]
}
krylov={
	'file':"kanamori_offdiag.krylov.h5",
	'name':"Krylov",
	'path':["G_%(spin)s"]
}

objects_to_plot = [[hyb,krylov],[hyb,krylov_qn],[ed,krylov],[ed,krylov_qn],[hyb,ed]]

num_o = 2

pp = PdfPages('G.pdf')
for plot_objs in objects_to_plot:
	try:
		plt.clf()
		for obj in plot_objs:
			arch = HDFArchive(obj['file'],'r')

			target_up = arch
			for p in obj['path']: target_up = target_up[p % {'spin':'up'}]    
			for o1,o2 in product(range(0,num_o),range(0,num_o)):
			    oplot(target_up[o1,o2], name=obj['name'] + (",$\uparrow%i%i$"%(o1,o2)))
			
			target_down = arch
			for p in obj['path']: target_down = target_down[p % {'spin':'down'}]
			for o1,o2 in product(range(0,num_o),range(0,num_o)):
			    oplot(target_down[o1,o2], name=obj['name'] + (",$\downarrow%i%i$"%(o1,o2)))
			
		setup_fig()
		pp.savefig(plt.gcf())

	except IOError: pass
	
pp.close()


