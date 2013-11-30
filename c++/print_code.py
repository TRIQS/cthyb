sources = """
. operator.hpp fock_state.hpp state.hpp fundamental_operator_set.hpp complete_hilbert_space.hpp partial_hilbert_space.hpp imperative_operator.hpp sorted_spaces.hpp krylov_worker.hpp exp_h_worker.hpp configuration.hpp qmc_data.hpp atomic_correlators_worker.hpp move_change_boundary_state.hpp move_insert.hpp move_remove.hpp measure_boundary_state.hpp measure_g.hpp measure_z.hpp statistics.hpp ctqmc_krylov.cpp ctqmc_krylov.hpp
"""

src_dir = "/home/parcolle/triqs/src/cthyb_krylov/c++"
code_filename = "~/ctkrylov_code.ps"

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


