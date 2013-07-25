@preamble@
from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages
import itertools

def setup_fig():
	axes = plt.gca()
	axes.set_ylabel('$G(\\tau)$')
	axes.legend(loc='lower center',prop={'size':10})

hyb={	
	'file':"spinless.hyb.h5",
     	'name':"HYB",
     	'path':["G_tau","tot"]
}
ed={	
	'file':"spinless.ed.h5",
	'name':"ED",
	'path':["G_tau"]
}
krylov_qn={
	'file':"spinless.krylov_qn.h5",
	'name':"Krylov (QN)",
	'path':["G_tau"],
}
krylov={
	'file':"spinless.krylov.h5",
	'name':"Krylov",
	'path':["G_tau"],
}

objects_to_plot=	[[hyb,krylov],[hyb,krylov_qn],
			[ed,krylov],[ed,krylov_qn],
			[hyb,ed]]

pp = PdfPages('G.pdf')
for plot_objs in objects_to_plot:
	try:
		plt.clf()
		for obj in plot_objs:
			arch = HDFArchive(obj['file'],'r')
			
			target = arch
			for p in obj['path']: target = target[p]
			
			for indices in itertools.product((0,1),(0,1)):
				oplot(target[indices[0],indices[1]], name=obj['name'] + ",%i%i" % indices)
			
		setup_fig()
		pp.savefig(plt.gcf())

	except IOError: pass

pp.close()
