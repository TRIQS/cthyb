#!/bin/env python

stats_filename = "krylov.stats.dat"

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

stats = {}
for line in open(stats_filename,'r'):
	if line[0] == '#': continue
	space_dim, krylov_dim, count = map(lambda s: int(s), line.split('\t'))
	stats[(space_dim,krylov_dim)] = count

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x=[]
y=[]
z=[]

for nm in stats:
    x.append(nm[0])
    y.append(nm[1])
    z.append(stats[nm])

dx=0.9
dy=0.9
dz=1.0

ax.bar3d(x,y,z,dx,dy,dz)
ax.set_xlabel("Space Dim")
ax.set_ylabel("Krylov Dim")
ax.set_zlabel("Counts")

plt.show()
