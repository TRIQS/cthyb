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
Utility functions
"""
from math import ceil
from numpy import argmax

def estimate_nfft_buf_size(gf_struct, pert_order_histograms):
    buf_sizes = {}
    for bn in gf_struct:
        if not bn in pert_order_histograms:
            raise RuntimeError("estimate_nfft_buf_size: no histogram for block '%s' is provided" % bn)
        else:
            max_order = argmax(pert_order_histograms[bn].data)
            block_size = len(gf_struct[bn])
            buf_sizes[bn] = int(max(ceil((max_order * max_order) / (block_size * block_size)), 1))

    return buf_sizes
