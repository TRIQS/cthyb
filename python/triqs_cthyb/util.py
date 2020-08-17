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

def multiplet_analysis(rho, h_loc_diag, orb_names, spin_names=['up','down'], off_diag=True):
    r"""
    Computes operator expectation values from measured
    density matrix and h_loc_diag object from cthyb solver.

    Measures N, Sz, S(S+1), constructs eigenstates in Fock state
    basis, assigns the according probabilities from rho for all
    impurity eigenstates, and stores the data in a Panda DataFrame.
    For more information check the guide in the documentation of cthyb
    regarding `Multiplet analysis & particle number histograms`.

    Warning: function assumes that N, Sz, and S(S+1) are good
    quantum numbers.

    Parameters
    ----------
    rho : numpy array
        measured density matrix from cthyb or equivalent solver, structured in subspaces
    h_loc_diag: triqs atom_diag object
        contains information about the local Hamiltonian (Hloc_0 + H_int)
    orb_names : list of int
        list of orbital indices
    spin_names : list of string
        list of strings containing the spin channel names
    off_diag: boolean
        determines whether blocks of Gf are named up_0 (false) or just up (true)

    Returns:
    --------
    res : Panda DataFrame
        containing all results structured
    """
    import pandas as pd

    from triqs.operators.util import make_operator_real
    from triqs.operators.util.observables import S_op, S2_op
    from triqs.atom_diag import quantum_number_eigenvalues
    from triqs.operators import n

    # res will be a list of dictionaries to be a panda data frame
    res = []

    # get fundamental operators from atom_diag object
    occ_operators = [n(*op) for op in h_loc_diag.fops]

    # construct total occupation operator from list
    N_op = sum(occ_operators)

    # create S2 and Sz operator
    S2 = S2_op(spin_names, orb_names, off_diag=off_diag)
    S2 = make_operator_real(S2)

    Sz=S_op('z', spin_names, orb_names, off_diag=off_diag)
    Sz = make_operator_real(Sz)

    # get eigenvalues
    S2_states = quantum_number_eigenvalues(S2, h_loc_diag)
    Sz_states = quantum_number_eigenvalues(Sz, h_loc_diag)

    # get particle numbers from h_loc_diag
    particle_numbers = quantum_number_eigenvalues(N_op, h_loc_diag)

    N_max = int(max(map(max, particle_numbers)))

    for sub in range(0,h_loc_diag.n_subspaces):

        # first get Fock space spanning the subspace
        fs_states = []
        for ind, fs in enumerate(h_loc_diag.fock_states[sub]):
            # get state in binary representation
            state = bin(int(fs))[2:].rjust(N_max, '0')
            fs_states.append("|"+state+">")

        # now extract particle number etc.
        for ind in range(h_loc_diag.get_subspace_dim(sub)):

            # get particle number
            # carefully here to not cast to int as the particle number
            # can be something like 1.999999996 and would get then 1!
            particle_number = round(particle_numbers[sub][ind])
            if abs(particle_number-particle_numbers[sub][ind]) > 1e-8:
                raise ValueError('round error for particle number to large!',
                                 particle_numbers[sub][ind])
            else:
                particle_number = int(particle_number)
            # energy of state
            eng=h_loc_diag.energies[sub][ind]

            # construct eigenvector in Fock state basis:
            ev_state = ''
            for i, elem in enumerate(h_loc_diag.unitary_matrices[sub][:,ind]):
                ev_state += ' {:+1.4f}'.format(elem)+fs_states[i]

            # get spin state
            ms=Sz_states[sub][ind]
            s_square=S2_states[sub][ind]

            # add to dict which becomes later the pandas data frame
            res.append({"Sub#" : sub,
                        "EV#" : ind,
                        "N" : particle_number,
                        "energy" : eng,
                        "prob": rho[sub][ind,ind],
                        "S2": abs(round(s_square,2)),
                        "m_s": round(ms,1),
                        "|m_s|": abs(round(ms,1)),
                        "state": ev_state})
    # panda data frame from res
    res = pd.DataFrame(res, columns=res[0].keys())

    return res
