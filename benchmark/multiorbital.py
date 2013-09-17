#!/bin/env python

from math import sqrt
from scipy.misc import factorial as fact
from itertools import product
import numpy as np

# Wigner 3-j symbols
def three_j_symbol(jm1, jm2, jm3):
    j1, m1 = jm1
    j2, m2 = jm2
    j3, m3 = jm3
    
    if (m1+m2+m3 != 0 or
        m1 < -j1 or m1 > j1 or
        m2 < -j2 or m2 > j2 or
        m3 < -j3 or m3 > j3 or
        j3 > j1 + j2 or
        j3 < abs(j1-j2)):
        return .0
        
    result = -1.0 if (j1-j2-m3) % 2 else 1.0
    result *= sqrt(fact(j1+j2-j3)*fact(j1-j2+j3)*fact(-j1+j2+j3)/fact(j1+j2+j3+1))
    result *= sqrt(fact(j1-m1)*fact(j1+m1)*fact(j2-m2)*fact(j2+m2)*fact(j3-m3)*fact(j3+m3))

    t_min = max(j2-j3-m1,j1-j3+m2,0)
    t_max = min(j1-m1,j2+m2,j1+j2-j3)
    
    t_sum = 0
    for t in range(t_min,t_max+1):
        t_sum += (-1.0 if t % 2 else 1.0)/(fact(t)*fact(j3-j2+m1+t)*fact(j3-j1-m2+t)*fact(j1+j2-j3-t)*fact(j1-m1-t)*fact(j2+m2-t))
    
    result *= t_sum
    return result

# Clebsch-Gordan coefficients
def clebsch_gordan(jm1, jm2, jm3):
    norm = sqrt(2*jm3[0]+1)*(-1 if jm1[0]-jm2[0]+jm3[1] % 2 else 1)
    return norm*three_j_symbol(jm1,jm2,(jm3[0],-jm3[1]))

# Angular matrix elements of particle-particle interaction
def angular_matrix_element(l, k, mp1, mp2, m2, m1):
    result = 0
    for q in range(-k,k+1):
        result += three_j_symbol((l,-mp1),(k,-q),(l,m2))*three_j_symbol((l,-mp2),(k,q),(l,m1))*(-1.0 if (mp1+q+mp2) % 2 else 1.0)
    result *= (2*l+1)**2 * (three_j_symbol((l,0),(k,0),(l,0))**2)
    return result

# The interaction matrix in the basis of spherical Harmonics
def U_matrix_spherical(radial_integrals, slater_notation=False):
    L = len(radial_integrals)-1
    if L<0: raise ValueError("radial_integrals list must not be empty.")

    # Full interaction matrix
    # Basis of spherical harmonics Y_{-2}, Y_{-1}, Y_{0}, Y_{1}, Y_{2}
    U_matrix = np.zeros((2*L+1,2*L+1,2*L+1,2*L+1), dtype=float)

    m_range = range(-L,L+1)
    for n, F in enumerate(radial_integrals):
        k = 2*n
        for mp1, mp2, m1, m2 in product(m_range,m_range,m_range,m_range):
            if slater_notation:
                # Slater notation
                U_matrix[mp1+L,mp2+L,m2+L,m1+L] += F * angular_matrix_element(L,k,mp1,mp2,m2,m1)
            else:
                # m1 <-> m2: 'Field theory' notation
                U_matrix[mp1+L,mp2+L,m1+L,m2+L] += F * angular_matrix_element(L,k,mp1,mp2,m2,m1)

    return U_matrix

# Transform the interaction matrix into another basis
def transform_U_matrix(U_matrix, W):
    # Transform the U-matrix
    return np.einsum("ij,kl,jlmo,mn,op",np.conj(W),np.conj(W),U_matrix,np.transpose(W),np.transpose(W))

def spherical2cubic(L):
    if L == 0:
        cubic_names = ("s",)
        W = np.array([1],dtype=complex)
    elif L == 1:
        cubic_names = ("x","y","z")
        W = np.zeros((3,3),dtype=complex)
        W[0,0] = 1.0/sqrt(2);   W[0,2] = -1.0/sqrt(2)
        W[1,0] = 1j/sqrt(2.0);  W[1,2] = 1j/sqrt(2.0)
        W[2,2] = 1.0
    elif L == 2:
        cubic_names = ("xy","yz","z^2","xz","x^2-y^2")
        W = np.zeros((5,5),dtype=complex)
        W[0,0] = 1j/sqrt(2);    W[0,4] = - 1j/sqrt(2)
        W[1,1] = 1j/sqrt(2);    W[1,3] = 1j/sqrt(2)
        W[2,2] = 1.0
        W[3,1] = 1.0/sqrt(2);   W[3,3] = -1.0/sqrt(2)
        W[4,0] = 1.0/sqrt(2);   W[4,4] = 1.0/sqrt(2)
    else: raise ValueError("spherical2cubic: implemented only for L=0,1,2")

    return (cubic_names,W)

spin_names = ("up","dn")
