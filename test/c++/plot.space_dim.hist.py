#!/bin/env python

stats_filename = "krylov.stats.dat"

import matplotlib.pyplot as plt
import numpy as np

stats = {}

for line in open(stats_filename,'r'):
    if line[0] == '#': continue
    space_dim, krylov_dim, count = map(lambda s: int(s), line.split('\t'))
    if space_dim in stats:
        stats[space_dim] += count
    else:
        stats[space_dim] = count

dims = sorted(stats)
counts = [stats[d] for d in dims]

ind = np.arange(len(stats))
width = 0.5

plt.title('Subspace dimension statistics')
plt.xlabel('Subspace dimension')
plt.ylabel('Times $\\exp(-\\hat H\\tau)$ computed')
plt.xticks(ind+width/2., dims)

plt.bar(ind,counts)

plt.show()
