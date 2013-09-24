#!/bin/env python

stats_filename = "dims.stats.dat"

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

stats = {}
space_dims = []
total_counts = 0
average_krylov_dim = 0
for line in open(stats_filename,'r'):
	if line[0] == '#': continue
	space_dim, krylov_dim, count = map(lambda s: int(s), line.split('\t'))
	
	if not space_dim in space_dims:
		space_dims.append(space_dim)
		stats[space_dim] = {}
	stats[space_dim][krylov_dim] = count
	
	average_krylov_dim += krylov_dim*count
	total_counts += count

space_dims.sort()

fig = plt.gcf()
ax = fig.add_subplot(111, projection='3d')

for n, space_dim in enumerate(space_dims):
	counts = np.zeros(space_dim,dtype=int)
	s = stats[space_dim]
	for krylov_dim in s: counts[krylov_dim-1] = s[krylov_dim]
	ax.bar(np.arange(0,space_dim),counts,zs=0.1*n,zdir='x')

ax.set_xticklabels(space_dims)
ax.set_xlabel("Space Dim")
ax.set_ylabel("Krylov Dim")
ax.set_zlabel("Count")

print "Average Krylov dimension:", float(average_krylov_dim)/total_counts

plt.show()
