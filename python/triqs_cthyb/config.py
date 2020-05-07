import h5py
import numpy as np
from triqs.plot.mpl_interface import *

# For the moment this is written for a single band case -- can be easily generalised.
class configuration:

    def __init__(self, oplist):
        self.oplist = oplist

    def plot(self, beta, conf_offset=0.0):
        fill = lambda dagger : 'full' if dagger else 'none'
        shape = lambda block : '<' if block else 'o'
        color = lambda block : 'r' if block else 'b'
        block_offset = lambda block : 0.003 if block else 0.0
        for i in range(self.oplist.shape[0]):
            plt.plot(self.oplist[i,0], [conf_offset + block_offset(self.oplist[i,1])], marker=shape(self.oplist[i,1]), color=color(self.oplist[i,2]), fillstyle=fill(self.oplist[i,3]))
            plt.hlines(conf_offset,beta,0.0)
        plt.xlim(beta,0.0)

def load_configuration(hdf_file,cid):
    hgroup = hdf_file['c_'+str(cid)]
    oplist = []
    for tau in hgroup: # list of taus in the given config
        block = hgroup[tau]['block'].value
        inner = hgroup[tau]['inner'].value
        dagger = hgroup[tau]['dagger'].value
        oplist.append([float(tau),block,inner,dagger])
    oplist = np.array(oplist)
    return configuration(oplist)

def count_configs(hdf_file,n_configs):
    empty = 0; up = 0; down = 0 ; updown = 0
    for i in range(1,n_configs+1):
        conf = load_configuration(hdf_file,i)
        if len(conf.oplist) == 0:
            empty += 1
        elif set(conf.oplist[:,1]) == set([0.0]):
            up += 1
        elif set(conf.oplist[:,1]) == set([1.0]):
            down += 1
        elif set(conf.oplist[:,1]) == set([0.0,1.0]):
            updown += 1
    print('empty = ',empty)
    print('up = ',up)
    print('down = ',down)
    print('up and down = ',updown)

def plot_configs(hdf_file,beta,n_configs,delta_configs):
    conf_offset = 0.0
    for conf_idx in range(delta_configs,n_configs,delta_configs):
        conf = load_configuration(hdf_file,conf_idx)
        conf.plot(beta,conf_offset)
        conf_offset += 0.01
        
def hist_pert_order(hdf_file,n_configs):
    length=[]
    for i in range(1,n_configs+1):
        conf = load_configuration(hdf_file,i)
        length.append(len(conf.oplist))
    plt.hist(length,bins=list(range(0,20)))
