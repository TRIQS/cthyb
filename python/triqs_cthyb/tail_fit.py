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

import numpy as np

from triqs.gf.gf_fnt import fit_hermitian_tail_on_window, replace_by_tail

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
