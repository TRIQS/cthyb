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

from triqs.gf.gf_fnt import fit_hermitian_tail_on_window, replace_by_tail

from triqs.operators import c, c_dag
from triqs.atom_diag import trace_rho_op

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



def tail_fit(
        Sigma_iw,
        fit_min_n=None, fit_max_n=None,
        fit_min_w=None, fit_max_w=None,
        fit_max_moment=None, fit_known_moments=None
        ):
    """
    Fit a high frequency 1/(iw)^n expansion of Sigma_iw 
    and replace the high frequency part with the fitted high frequency expansion.

    Either give frequency window to fit on in terms of matsubara frequencies index 
    (fit_min_n/fit_max_n) or value (fit_min_w/fit_max_w).

    Parameters
    ----------
    Sigma_iw : Gf
               Self-energy.
    fit_min_n : int, optional, default=int(0.8*len(Sigma_iw.mesh))
                Matsubara frequency index from which tail fitting should start.
    fit_max_n : int, optional, default=int(len(Sigma_iw.mesh))
                Matsubara frequency index at which tail fitting should end.
    fit_min_w : float, optional
                Matsubara frequency from which tail fitting should start.
    fit_max_w : float, optional
                Matsubara frequency at which tail fitting should end.
    fit_max_moment : int, optional
                     Highest moment to fit in the tail of Sigma_iw.
    fit_known_moments : ``ndarray.shape[order, Sigma_iw[0].target_shape]``, optional, default = None
                        Known moments of Sigma_iw, given as an numpy ndarray

    Returns
    -------
    Sigma_iw : Gf
               Self-energy.
    """

    # Define default tail quantities
    if fit_min_w is not None: fit_min_n = int(0.5*(fit_min_w*Sigma_iw.mesh.beta/np.pi - 1.0))
    if fit_max_w is not None: fit_max_n = int(0.5*(fit_max_w*Sigma_iw.mesh.beta/np.pi - 1.0))
    if fit_min_n is None: fit_min_n = int(0.8*len(Sigma_iw.mesh)/2)
    if fit_max_n is None: fit_max_n = int(len(Sigma_iw.mesh)/2)
    if fit_max_moment is None: fit_max_moment = 3


    if fit_known_moments is None:
        fit_known_moments = {}
        for name, sig in Sigma_iw:
            shape = [0] + list(sig.target_shape)
            fit_known_moments[name] = np.zeros(shape, dtype=np.complex) # no known moments

    # Now fit the tails of Sigma_iw and replace the high frequency part with the tail expansion
    for name, sig in Sigma_iw:

        tail, err = fit_hermitian_tail_on_window(
            sig,
            n_min = fit_min_n,
            n_max = fit_max_n,
            known_moments = fit_known_moments[name],
            # set max number of pts used in fit larger than mesh size, to use all data in fit
            n_tail_max = 10 * len(sig.mesh), 
            expansion_order = fit_max_moment
            )
        
        replace_by_tail(sig, tail, n_min=fit_min_n)        

    return Sigma_iw
