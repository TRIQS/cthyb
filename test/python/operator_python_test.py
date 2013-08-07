#!/bin/env python

from pytriqs.applications.impurity_solvers.cthyb_krylov import *
import itertools

C_list = [C(1,0),C(2,0),C(3,0)]
Cd_list = [C_dag(1,0), C_dag(2,0), C_dag(3,0)]

print "Anticommutators:"
for cd,c in itertools.product(Cd_list,C_list):
	print "{", cd, ",", c, "} =", cd*c + c*cd

print "Commutators:"
for cd,c in itertools.product(Cd_list,C_list): 
	print "[", cd, ",", c, "] =", cd*c - c*cd


# Algebra
x = C(0,0)
y = C_dag(1,0)

print
print "Algebra:"  
print "x =", x
print "y =", y

print "-x =", -x
print "x + 2.0 =", x + 2.0
print "2.0 + x =", 2.0 + x
print "x - 2.0 =", x - 2.0
print "2.0 - x =", 2.0 - x
print "3.0*y =", 3.0*y
print "y*3.0 =", y*3.0
print "x + y =", x + y
print "x - y =", x - y
print "(x + y)*(x - y) =", (x + y)*(x - y)

# N^3
print
print "N^3:"
N = N(0,'up') + N(0,'dn')
N3 = N*N*N
print "N =", N
print "N^3 =", N3
