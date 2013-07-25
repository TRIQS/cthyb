@preamble@
from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import plt, oplot
from matplotlib.backends.backend_pdf import PdfPages

def setup_fig():
	axes = plt.gca()
	axes.set_ylabel('$G(\\tau)$')
	axes.legend(loc='lower center',prop={'size':10})

hyb={	
	'file':"anderson.hyb.h5",
     	'name':"HYB",
     	'up_path':["G_tau","up"],
     	'down_path':["G_tau","down"]
}
ed={	
	'file':"anderson.ed.h5",
	'name':"ED",
	'up_path':["G_up"],
	'down_path':["G_down"]
}
krylov_blocks_qn={
	'file':"anderson.krylov_block_qn.h5",
	'name':"Krylov (block,QN)",
	'up_path':["G_up"],
	'down_path':["G_down"]
}
krylov_blocks={
	'file':"anderson.krylov_block.h5",
	'name':"Krylov (block)",
	'up_path':["G_up"],
	'down_path':["G_down"]
}
krylov_qn={
	'file':"anderson.krylov_qn.h5",
	'name':"Krylov (QN)",
	'up_path':["G_tot",(0,0)],
	'down_path':["G_tot",(1,1)]
}
krylov={
	'file':"anderson.krylov.h5",
	'name':"Krylov",
	'up_path':["G_tot",(0,0)],
	'down_path':["G_tot",(1,1)]
}

objects_to_plot=	[[hyb,krylov],[hyb,krylov_qn],[hyb,krylov_blocks],[hyb,krylov_blocks_qn],
			[ed,krylov],[ed,krylov_qn],[ed,krylov_blocks],[ed,krylov_blocks_qn],
			[hyb,ed]]

pp = PdfPages('G.pdf')
for plot_objs in objects_to_plot:
	try:
		plt.clf()
		for obj in plot_objs:
			arch = HDFArchive(obj['file'],'r')

			target_up = arch
			for p in obj['up_path']: target_up = target_up[p]
			oplot(target_up, name=obj['name'] + ",$\uparrow\uparrow$")
			
			target_down = arch
			for p in obj['down_path']: target_down = target_down[p]
			oplot(target_down, name=obj['name'] + ",$\downarrow\downarrow$")
			
		setup_fig()
		pp.savefig(plt.gcf())

	except IOError: pass

pp.close()
