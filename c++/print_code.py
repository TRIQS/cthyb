sources = """
. qmc_parameters.hpp t_ordered_colored_c_ops.hpp segment.hpp trace_c_ops.hpp trace_c_ops.cpp configuration.hpp hybridization_dets.hpp hyb_opmap_no_storage.hpp move_insert_segment.hpp move_remove_segment.hpp move_move_segment.hpp move_swap_empty_lines.hpp ctqmc_seg.hpp ctqmc_seg.cpp measure_hist.hpp measure_gt.hpp measure_nn.hpp 
"""

src_dir = "~/triqs/src/applications/impurity_solvers/ctqmc_seg/"
code_filename = "~/ctseg_code.ps"

#################################################

import sys,os,tempfile
rep = tempfile.mkdtemp()
old_rep = os.getcwd()
os.chdir(rep)

print rep

def cut_license(f) : 
    status = 0
    print "treating %s"%f
    A = open(f).readlines()
    out = open (f,'w')
    for line in A :
      if line [0:3] == "/**" : status = 1
      if line.rstrip() [-3:] == "**/" : status = 3
      if status ==2 : out.write( line )
      if status ==3 : 
          status = 2

flist = []
for lines in sources.split('\n') : 
    sp = lines.split()
    if sp : 
        flist += [ sp[0] + '/' + x for x in sp[1:]]

print flist

flist2 = ' '.join(flist)
os.system (" cd %s && tar cvf %s/all.tar %s && cd %s && tar xvf all.tar"%(src_dir,rep,flist2,rep))

for f in flist : 
    cut_license(f)

os.system( "enscript -G -f Courier7 --color -Ecpp -o %s %s"%(code_filename,flist2)  )
os.chdir(old_rep)
os.system ("rm -rf %s"%rep)  


