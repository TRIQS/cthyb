################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2014 by P. Seth, I. Krivenko, M. Ferrero, O. Parcollet
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
r"""
CTHYB utility functions:
"""
from math import ceil
from numpy import argmax

import numpy as np
from triqs.operators import c, c_dag
from triqs.atom_diag import trace_rho_op

def block_size_from_gf_struct(block_name, gf_struct):
    bns, idxs = list(zip(*gf_struct))
    bidx = bns.index(block_name)
    block_size = len(idxs[bidx])
    return block_size

def estimate_nfft_buf_size(gf_struct, pert_order_histograms):
    buf_sizes = {}
    for bn, idxs in gf_struct:
        if not bn in pert_order_histograms:
            raise RuntimeError("estimate_nfft_buf_size: no histogram for block '%s' is provided" % bn)
        else:
            max_order = argmax(pert_order_histograms[bn].data)
            block_size = block_size_from_gf_struct(bn, gf_struct)
            buf_sizes[bn] = int(max(ceil((max_order * max_order) / (block_size * block_size)), 1))

    return buf_sizes

def orbital_occupations(density_matrix, gf_struct, h_loc_diag):
    
    dtype=density_matrix[0].dtype
    occ_mat = {bl: np.zeros((bl_size,bl_size), dtype=dtype) for bl, bl_size in gf_struct}

    for bl, bl_size in gf_struct:
        for i in range(bl_size):
            for j in range(bl_size):
                occ_mat[bl][i,j] = trace_rho_op(density_matrix, c_dag(bl, i)*c(bl,j), h_loc_diag)

    return occ_mat

