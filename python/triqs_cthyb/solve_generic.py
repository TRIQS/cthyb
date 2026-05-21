# Copyright (c) 2024-2026 Simons Foundation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You may obtain a copy of the License at
#     https://www.gnu.org/licenses/gpl-3.0.txt
#
# Authors: Nils Wentzell

"""
CT-HYB solver: generic functional interface with optional post-processing.

Provides a functional API to the triqs_cthyb solver with:
- Dynamic dispatch based on mesh type (MeshImFreq, MeshDLRImFreq)
- Optional post-processing strategies for self-energy extraction
- Always-on density matrix measurement for moment access
"""

from __future__ import annotations
from dataclasses import dataclass

import numpy as np
import triqs.utility.mpi as mpi
from triqs.gfs import (
    Gf, BlockGf, MeshImFreq, MeshDLRImFreq,
    make_gf_from_fourier, fit_hermitian_tail, make_hermitian,
    make_gf_imfreq, make_gf_dlr_imfreq, make_gf_imtime, make_gf_dlr_imtime,
    make_gf_dlr, fit_gf_dlr, inverse, iOmega_n,
)
from triqs.gfs.tools import make_zero_tail
from triqs.gfs.dlr_crm_dyson_solver import minimize_dyson
from triqs.operators import Operator
from triqs.operators.util import op_from_block_matrix, block_matrix_from_op

from triqs.solver_utils import SolverResults, tail_fit
from triqs_cthyb.solver import Solver


# =============================================================================
# Post-Processing Strategy Parameter Dataclasses
# =============================================================================

@dataclass(frozen=True)
class TailFitParams:
    """Parameters for high-frequency tail fitting.

    Fits a 1/(iw)^n expansion to the high-frequency part of Sigma(iw) and
    replaces the noisy high-frequency data with the fitted tail.

    Either specify the fitting window via Matsubara indices (fit_min_n, fit_max_n)
    or via frequency values (fit_min_w, fit_max_w). Frequency values take precedence.
    """
    fit_min_n: int | None = None
    fit_max_n: int | None = None
    fit_min_w: float | None = None
    fit_max_w: float | None = None
    fit_max_moment: int = 3


@dataclass(frozen=True)
class LegendreParams:
    """Parameters for Legendre polynomial filtering/fitting.

    Can use either measured Legendre coefficients from solver (use_measured=True)
    or fit G(tau) to Legendre polynomials after the solve (use_measured=False).
    """
    n_l: int = 40
    use_measured: bool = False


@dataclass(frozen=True)
class CRMParams:
    """Parameters for Constrained Residual Minimization (CRM) Dyson solver.

    Uses DLR to fit G(tau), then minimizes the Dyson equation residual to extract Sigma.
    """
    w_max: float = 10.0
    eps: float = 1e-10


PostProcessParams = TailFitParams | LegendreParams | CRMParams


# =============================================================================
# Legendre Filter
# =============================================================================

def _legendre_filter(G_tau: BlockGf, n_l: int) -> BlockGf:
    """Fit G(tau) to Legendre polynomials to filter high-frequency noise."""
    from triqs.gfs import MeshLegendre, GfLegendre

    beta = G_tau.mesh.beta
    names, blocks = [], []
    for block, g_tau in G_tau:
        mesh_l = MeshLegendre(beta=beta, S='Fermion', n_max=n_l)
        g_l = GfLegendre(mesh=mesh_l, target_shape=g_tau.target_shape)
        g_l.set_from_imtime(g_tau)
        names.append(block)
        blocks.append(g_l)
    return BlockGf(name_list=names, block_list=blocks, make_copies=False)


def _set_gf_from_legendre(G_l, G_iw, G_tau):
    """Set G(iw) and G(tau) from Legendre representation (in-place)."""
    G_l.enforce_discontinuity(np.identity(G_l.target_shape[0]))
    G_iw.set_from_legendre(G_l)
    G_tau.set_from_legendre(G_l)


# =============================================================================
# Post-Processing Strategy Functions
# =============================================================================

def _dyson_result(solver, G_iw, G_tau, Sigma_iw):
    """Build the standard post-processing result dict from Sigma_iw."""
    Sigma_dynamic = Sigma_iw.copy()
    for bl, g in Sigma_dynamic:
        g << g - solver.Sigma_Hartree[bl]

    return {
        'G_iw': G_iw,
        'G_tau': G_tau,
        'Sigma_iw': Sigma_iw,
        'Sigma_dynamic': Sigma_dynamic,
        'Sigma_HartreeFock': [solver.Sigma_Hartree[bl] for bl, _ in solver.G_iw],
    }


def _postprocess_dyson(solver, G0_iw):
    """Extract self-energy via plain Dyson equation: Sigma = G0^{-1} - G^{-1}."""
    G_iw = make_hermitian(solver.G_iw.copy())
    G_tau = solver.G_tau.copy()

    Sigma_iw = G0_iw.copy()
    Sigma_iw << inverse(G0_iw) - inverse(G_iw)

    return _dyson_result(solver, G_iw, G_tau, Sigma_iw)


def _postprocess_tail_fit(solver, G0_iw, params):
    """Extract self-energy via Dyson equation followed by tail fitting."""
    G_iw = make_hermitian(solver.G_iw.copy())
    G_tau = solver.G_tau.copy()

    Sigma_iw = G0_iw.copy()
    Sigma_iw << inverse(G0_iw) - inverse(G_iw)

    known_moments = dict(solver.Sigma_moments) if solver.Sigma_moments else None
    tail_fit(
        Sigma_iw,
        fit_min_n=params.fit_min_n,
        fit_max_n=params.fit_max_n,
        fit_min_w=params.fit_min_w,
        fit_max_w=params.fit_max_w,
        fit_max_moment=params.fit_max_moment,
        fit_known_moments=known_moments,
    )

    return _dyson_result(solver, G_iw, G_tau, Sigma_iw)


def _postprocess_legendre(solver, G0_iw, params):
    """Extract self-energy via Legendre filtering then Dyson equation."""
    if params.use_measured:
        if not hasattr(solver, 'G_l') or solver.G_l is None:
            raise ValueError("use_measured=True but solver has no G_l. "
                           "Set measure_G_l=True in solver params.")
        G_l = solver.G_l.copy()
    else:
        G_l = _legendre_filter(solver.G_tau, params.n_l)

    G_iw = G0_iw.copy()
    G_tau = solver.G_tau.copy()

    for bl, g_l in G_l:
        _set_gf_from_legendre(g_l, G_iw[bl], G_tau[bl])

    G_iw = make_hermitian(G_iw)

    Sigma_iw = G0_iw.copy()
    Sigma_iw << inverse(G0_iw) - inverse(G_iw)

    result = _dyson_result(solver, G_iw, G_tau, Sigma_iw)
    result['G_l'] = G_l
    return result


def _postprocess_crm(solver, G0_iw, params):
    """Extract self-energy via Constrained Residual Minimization (CRM)."""
    result = None

    if mpi.is_master_node():
        mpi.report(f'\nCRM Dyson solver with (w_max, eps) = ({params.w_max}, {params.eps})\n')

        G_dlr = fit_gf_dlr(solver.G_tau, w_max=params.w_max, eps=params.eps)
        G_tau_dlr = make_gf_dlr_imtime(G_dlr)
        G_iw_dlr = make_gf_dlr_imfreq(G_dlr)

        mesh_dlr_iw = MeshDLRImFreq(G_dlr.mesh)
        G0_dlr_iw = G0_iw.copy()
        for bl, gf in G0_iw:
            # Project G0 onto mesh_dlr_iw by going through DLR coefficients.
            # DLR coeffs can be evaluated at arbitrary Matsubara frequencies.
            g_dlr = make_gf_dlr(gf) if isinstance(gf.mesh, MeshDLRImFreq) else gf
            g_proj = Gf(mesh=mesh_dlr_iw, target_shape=gf.target_shape)
            for iwn in mesh_dlr_iw:
                g_proj[iwn] = g_dlr(iwn)
            G0_dlr_iw[bl] = g_proj

        Sigma_dlr = G0_dlr_iw.copy()
        np.random.seed(85281)

        for bl, gf in Sigma_dlr:
            mpi.report(f'Minimizing Dyson via CRM for Sigma[{bl}]')
            Sigma_dlr[bl], _, _ = minimize_dyson(
                G0_dlr=G0_dlr_iw[bl],
                G_dlr=G_iw_dlr[bl],
                Sigma_moments=solver.Sigma_moments[bl]
            )

        n_iw = len(G0_iw.mesh) // 2 if hasattr(G0_iw.mesh, '__len__') else solver.n_iw
        Sigma_iw = make_gf_imfreq(Sigma_dlr, n_iw=n_iw)

        for bl, gf in Sigma_iw:
            gf += solver.Sigma_moments[bl][0]

        Sigma_dynamic = Sigma_dlr.copy()

        G_iw = make_hermitian(solver.G_iw.copy())
        G_tau = solver.G_tau.copy()

        result = {
            'G_iw': G_iw,
            'G_tau': G_tau,
            'G_tau_dlr': G_tau_dlr,
            'Sigma_iw': Sigma_iw,
            'Sigma_dynamic': Sigma_dynamic,
            'Sigma_dlr': Sigma_dlr,
            'Sigma_HartreeFock': [solver.Sigma_Hartree[bl] for bl, _ in solver.G_iw],
        }

    mpi.barrier()
    result = mpi.bcast(result)

    return result


def _apply_postprocessing(solver, G0_iw, postprocess):
    """Apply post-processing strategy to extract self-energy from solver output."""
    if postprocess == 'dyson':
        return _postprocess_dyson(solver, G0_iw)
    if isinstance(postprocess, TailFitParams):
        return _postprocess_tail_fit(solver, G0_iw, postprocess)
    if isinstance(postprocess, LegendreParams):
        return _postprocess_legendre(solver, G0_iw, postprocess)
    if isinstance(postprocess, CRMParams):
        return _postprocess_crm(solver, G0_iw, postprocess)
    raise TypeError(f"Unknown postprocess type: {type(postprocess)}. "
                   f"Expected 'dyson' or one of: TailFitParams, LegendreParams, CRMParams")


# =============================================================================
# Internal Helpers
# =============================================================================

def _extract_interface_params(solver_params):
    """Extract and remove interface-specific params from solver_params."""
    params = solver_params.copy()
    n_tau = params.pop('n_tau', 10001)
    n_l = params.pop('n_l', 40)
    return n_tau, n_l, params


def _canonicalize_h_loc0(h_loc0, gf_struct):
    """Normalize h_loc0 into (Operator, block-matrix array) regardless of input form.

    Accepts either a many_body_operator or an iterable of dense block matrices
    (one per block in gf_struct).
    """
    if isinstance(h_loc0, Operator):
        return h_loc0, block_matrix_from_op(h_loc0, gf_struct)
    bl = np.empty(len(h_loc0), dtype=object)
    for i, m in enumerate(h_loc0):
        bl[i] = np.asarray(m)
    return op_from_block_matrix(bl, gf_struct), bl


def _prepare_solver(gf_struct, beta, n_iw, n_tau, n_l, solver_params, postprocess):
    """Prepare the CT-HYB solver instance and parameters."""
    solver_params = solver_params.copy()
    solver_params['measure_density_matrix'] = True
    solver_params['use_norm_as_weight'] = True

    if isinstance(postprocess, LegendreParams) and postprocess.use_measured:
        solver_params['measure_G_l'] = True

    S = Solver(
        gf_struct=gf_struct,
        beta=beta,
        n_iw=n_iw,
        n_tau=n_tau,
        n_l=n_l,
        delta_interface=True,
    )

    return S, solver_params


def _prepare_delta_tau_imfreq(Delta_iw, S):
    """Fourier transform Delta(iw) to Delta(tau) with hermitian tail fitting."""
    def _fourier_with_tail(G):
        tail = fit_hermitian_tail(G, make_zero_tail(G, 1))[0]
        return make_gf_from_fourier(G, S.Delta_tau[0].mesh, tail)

    for block, _ in S.Delta_tau:
        S.Delta_tau[block] << _fourier_with_tail(Delta_iw[block])


def _prepare_delta_tau_dlr(Delta_iw, S):
    """Convert Delta from DLR mesh to tau mesh."""
    S.Delta_tau << make_gf_imtime(Delta_iw, S.n_tau)


def _build_G0_iw(Delta_iw, hloc0_bl_mat):
    """Construct G0(iw) = [iw - h_loc0 - Delta(iw)]^{-1}."""
    G0_iw = Delta_iw.copy()
    for idx, (bl, g) in enumerate(G0_iw):
        G0_iw[bl] << inverse(iOmega_n - hloc0_bl_mat[idx] - Delta_iw[bl])
    return G0_iw


# =============================================================================
# Main Solve Function
# =============================================================================

def solve_generic(
    Delta_iw,
    h_loc0,
    h_int,
    postprocess=None,
    **solver_params,
):
    """Solve the quantum impurity problem using CT-HYB.

    Parameters
    ----------
    Delta_iw : BlockGf
        Hybridization function Delta(iw). Mesh type determines preprocessing:
        - MeshImFreq: Fourier transform with tail fitting
        - MeshDLRImFreq: Direct DLR to tau conversion
    h_loc0 : Operator | list[np.ndarray]
        Local non-interacting Hamiltonian. Accepted as either a many_body_operator
        or an iterable of dense block matrices (one per block in Delta_iw).
    h_int : Operator
        Interaction Hamiltonian.
    postprocess : str | PostProcessParams | None, optional
        Post-processing strategy for self-energy extraction:
        - None: Return raw solver output
        - 'dyson': Plain Dyson equation
        - TailFitParams(...): Dyson equation followed by tail fitting
        - LegendreParams(...): Legendre filtering then Dyson
        - CRMParams(...): Constrained Residual Minimization
    **solver_params
        Parameters passed to cthyb solver. Key parameters:
        - n_tau : int - Number of tau points (default: 10001)
        - n_l : int - Number of Legendre polynomials (default: 40)
        - length_cycle : int - MC cycle length
        - n_cycles : int - Number of MC cycles per MPI rank
        - n_warmup_cycles : int - Warmup cycles

    Returns
    -------
    SolverResults
        Container with solver output.
    """
    n_tau, n_l, solver_params = _extract_interface_params(solver_params)

    mesh = Delta_iw.mesh
    gf_struct = [(bl, gf.target_shape[0]) for (bl, gf) in Delta_iw]
    h_loc0_op, h_loc0_bl = _canonicalize_h_loc0(h_loc0, gf_struct)

    # Determine n_iw
    n_iw = len(mesh) // 2 if isinstance(mesh, MeshImFreq) else int(2*mesh.values()[-1].n)

    # Prepare solver
    S, solver_params = _prepare_solver(
        gf_struct, mesh.beta, n_iw, n_tau, n_l, solver_params, postprocess
    )

    # Prepare Delta(tau) based on mesh type
    if isinstance(mesh, MeshImFreq):
        _prepare_delta_tau_imfreq(Delta_iw, S)
    elif isinstance(mesh, MeshDLRImFreq):
        _prepare_delta_tau_dlr(Delta_iw, S)
    else:
        raise NotImplementedError(f"Unsupported mesh type: {type(mesh)}")

    # Solve
    mpi.report(f"Solving impurity problem with CT-HYB (postprocess={postprocess})")
    S.solve(h_loc0=h_loc0_op, h_int=h_int, **solver_params)

    # Post-process or return raw results
    if postprocess is None:
        return SolverResults(Solver=S)

    # Build G0 for post-processing
    G0_iw = _build_G0_iw(Delta_iw, h_loc0_bl)

    # For MeshDLRImFreq, we need G0 on ImFreq mesh for some strategies
    if isinstance(mesh, MeshDLRImFreq) and not isinstance(postprocess, CRMParams):
        G0_iw = make_gf_imfreq(G0_iw, n_iw=n_iw)

    # Apply post-processing
    pp_result = _apply_postprocessing(S, G0_iw, postprocess)

    result_kwargs = {k: v for k, v in pp_result.items() if v is not None}
    result_kwargs['Solver'] = S

    return SolverResults(**result_kwargs)
