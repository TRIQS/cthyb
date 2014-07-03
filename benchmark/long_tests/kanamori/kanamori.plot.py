@preamble@
from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import *
from matplotlib.backends.backend_pdf import PdfPages

def setup_fig():
	axes = plt.gca()
	axes.set_ylabel('$G(\\tau)$')
	axes.legend(loc='lower center',prop={'size':10})

hyb={	
	'file':"kanamori.hyb.h5",
     	'name':"HYB",
     	'path':["G_tau","%(spin)s-%(orbital)i"],
}
ed={	
	'file':"kanamori.ed.h5",
	'name':"ED",
	'path':["G_%(spin)s-%(orbital)i"],
}
krylov_qn={
	'file':"kanamori.krylov_qn.h5",
	'name':"Krylov (QN)",
	'path':["G_%(spin)s-%(orbital)i"]
}
krylov={
	'file':"kanamori.krylov.h5",
	'name':"Krylov",
	'path':["G_%(spin)s-%(orbital)i"]
}

objects_to_plot = [[hyb,krylov],[hyb,krylov_qn],[ed,krylov],[ed,krylov_qn],[hyb,ed]]

num_o = 2

pp = PdfPages('G.pdf')
for plot_objs in objects_to_plot:
	try:
		plt.clf()
		for obj in plot_objs:
			arch = HDFArchive(obj['file'],'r')
			
			for o in range(0,num_o):
			    target_up = arch
			    for p in obj['path']: target_up = target_up[p % {'spin':'up','orbital':o}]    
			    oplot(target_up, name=obj['name'] + (",$\uparrow%i$"%o))
			    target_down = arch
			    for p in obj['path']: target_down = target_down[p % {'spin':'down','orbital':o}]    
			    oplot(target_down, name=obj['name'] + (",$\downarrow%i$"%o))
						
		setup_fig()
		pp.savefig(plt.gcf())

	except IOError: pass

pp.close()

