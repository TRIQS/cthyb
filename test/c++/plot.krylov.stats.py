#!/bin/env python

stats_filename = "dims.stats.dat"

import matplotlib.pyplot as plt
import numpy as np

stats = {}

x = []
y = []
for line in open(stats_filename,'r'):
	if line[0] == '#': continue
	space_dim, krylov_dim, count = map(lambda s: int(s), line.split('\t'))
	stats[(space_dim,krylov_dim)] = count
	
	x.append(space_dim)
	y.append(krylov_dim)

x = np.arange(min(x),max(x)+2)
y = np.arange(min(y),max(y)+2)
z = np.zeros((len(x)-1,len(y)-1),dtype=int)

for nm in stats:
    xi = nm[0] - x[0]
    yi = nm[1] - y[0]
    z[xi,yi] = stats[nm]

plt.pcolor(x,y,z.T,cmap="gray")
#plt.imshow(z,extent=[x[0],x[-1],y[0],y[-1]])

ax = plt.gca()
ax.set_xlim(x[0],x[-1])
ax.set_ylim(y[0],y[-1])
ax.set_xlabel("Space Dim")
ax.set_ylabel("Krylov Dim")

plt.show()
