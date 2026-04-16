################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2017 by H. UR Strand, P. Seth, I. Krivenko,
#                       M. Ferrero, O. Parcollet
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
tail fitting and high frequency moments
"""
import numpy as np

from triqs.operators import c, c_dag
from triqs.atom_diag import trace_rho_op
from triqs.solver_utils import tail_fit  # noqa: F401 -- re-export for backwards compatibility

def _comm(A,B): return A*B - B*A
def _anticomm(A,B): return A*B + B*A

def sigma_high_frequency_moments(density_matrix,
                           ad_imp, 
                           gf_struct, 
                           h_int):
    """
    Calculate the first and second high frequency moment of Sigma_iw
    following Rev. Mod. Phys. 83, 349 (2011). They read
    (0) Sigma_0       = -<{[Hint,c],c+}> (Hartree shift)
    (1) Sigma_1       =  <{[Hint,[Hint,c]],c+}> - Sigma_0^2,
    where Hint is the interaction Hamiltonian

    Parameters
    ----------
    density_matrix : list, np.ndarray
                     measured density matrix from TRIQS/CTHYB.
    ad_imp         : AtomDiag
                     h_loc_diagonalization from TRIQS/CTHYB.
    gf_struct      : List of pairs (str,int)
                     Block structure of Green's function.
    h_int          : triqs.operators.Operator
                     interaction Hamiltonian

    Returns
    -------
    sigma_moments  : dict, np.ndarray
                     first and second moments in a dict with the
                     same block strucutre of the TRIQS Gf object.
    """


    sigma_moments = {bl : np.zeros((2, bl_size, bl_size),dtype=complex) for bl, bl_size in gf_struct}
    for bl, bl_size in gf_struct:
        for orb1 in range(bl_size):
            for orb2 in range(bl_size):

                # Sigma_HF term
                op_HF = -_anticomm(_comm(h_int, c(bl,orb1)), c_dag(bl,orb2))
                sigma_moments[bl][0,orb1,orb2] = trace_rho_op(density_matrix, op_HF, ad_imp)

                # Sigma_1/iwn term
                op_iw = _anticomm(_comm(h_int, _comm(h_int, c(bl,orb1))), c_dag(bl,orb2))
                sigma_moments[bl][1,orb1,orb2] = trace_rho_op(density_matrix, op_iw, ad_imp) - sigma_moments[bl][0,orb1,orb2]**2

    return sigma_moments


def green_high_frequency_moments(density_matrix,
                           ad_imp, 
                           gf_struct, 
                           h_imp):
    """
    Calculate the first and second high frequency moment of G_iw
    following Rev. Mod. Phys. 83, 349 (2011). They read
    (0) G_0           =    0
    (1) G_1           =  <{c,c+}>
    (2) G_2           = -<{[H,c],c+}>
    where H is the impurity Hamiltonian (H = impurity levels + Hint).

    Parameters
    ----------
    density_matrix : list, np.ndarray
                     measured density matrix from TRIQS/CTHYB.
    ad_imp         : AtomDiag
                     h_loc_diagonalization from TRIQS/CTHYB.
    gf_struct      : List of pairs (str,int)
                     Block structure of Green's function.
    h_imp          : triqs.operators.Operator
                     impurity Hamiltonian   

    Returns
    -------
    green_moments  : dict, np.ndarray
                     first and second moments in a dict with the
                     same block strucutre of the TRIQS Gf object.
    """

    green_moments = {bl : np.zeros((3, bl_size, bl_size),dtype=complex) for bl, bl_size in gf_struct}
    for bl, bl_size in gf_struct:
        # G_0/iwn = 1/iwn
        green_moments[bl][1] = np.eye(bl_size)
        for orb1 in range(bl_size):
            for orb2 in range(bl_size):
                # G_1/iwn**2 term
                op = -_anticomm(_comm(h_imp, c(bl,orb1)), c_dag(bl,orb2))
                green_moments[bl][2,orb1,orb2] = trace_rho_op(density_matrix, op, ad_imp)

    return green_moments
