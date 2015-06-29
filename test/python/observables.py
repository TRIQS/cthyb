"""Operators for commonly used observables"""

import numpy as np
from math import sqrt
from pytriqs.operators.operators import Operator, n, c_dag, c
from pytriqs.operators.util.op_struct import get_mkind
from pytriqs.operators.util.U_matrix import spherical_to_cubic
from itertools import product

pauli_matrix = {'x' : np.array([[0,1],[1,0]]),
                'y' : np.array([[0,-1j],[1j,0]]),
                'z' : np.array([[1,0],[0,-1]]),
                '+' : np.array([[0,2],[0,0]]),
                '-' : np.array([[0,0],[2,0]])}

def N_op(spin_names, orb_names, off_diag = None, map_operator_structure = None):
    r"""
    Create an operator of the total number of particles.

    .. math:: \hat N = \sum_{i\sigma} a_{i\sigma}^\dagger a_{i\sigma}.

    Parameters
    ----------
    spin_names : list of strings
                 Names of the spins, e.g. ['up','down'].
    orb_names : list of strings or int
                Names of the orbitals, e.g. [0,1,2] or ['t2g','eg'].
    off_diag : boolean
               Do we have (orbital) off-diagonal elements?
               If yes, the operators and blocks are denoted by ('spin', 'orbital'),
               otherwise by ('spin_orbital',0).
    map_operator_structure : dict 
                             Mapping of names of GF blocks names from one convention to another, 
                             e.g. {('up', 0): ('up_0', 0), ('down', 0): ('down_0',0)}.
                             If provided, the operators and blocks are denoted by the mapping of ``('spin', 'orbital')``.

    Returns
    -------
    N : Operator
        The total number of particles.

    """
    mkind = get_mkind(off_diag,map_operator_structure)
    N = Operator()
    for sn, on in product(spin_names,orb_names): N += n(*mkind(sn,on))
    return N

def S_op(component, spin_names, orb_names, off_diag = None, map_operator_structure = None):
    r"""
    Create a component of the spin vector operator.

    .. math:: \hat S^{x,y,z} = \frac{1}{2}\sum_{i\sigma\sigma'}

    Parameters
    ----------
    component : string
                Component to be created, one of 'x','y','z','+', or '-'.
                WARNING: y-component is not supported at the moment!
    spin_names : list of strings
                 Names of the spins, e.g. ['up','down'].
    orb_names : list of strings or int
                Names of the orbitals, e.g. [0,1,2] or ['t2g','eg'].
    off_diag : boolean
               Do we have (orbital) off-diagonal elements?
               If yes, the operators and blocks are denoted by ('spin', 'orbital'),
               otherwise by ('spin_orbital',0).
    map_operator_structure : dict
                             Mapping of names of GF blocks names from one convention to another,
                             e.g. {('up', 0): ('up_0', 0), ('down', 0): ('down_0',0)}.
                             If provided, the operators and blocks are denoted by the mapping of ``('spin', 'orbital')``.

    Returns
    -------
    S : Operator
        The component of the spin vector operator.

    """
    # FIXME
    assert component != 'y', "We cannot construct operators with complex coefficients at the moment. Sorry for that!"

    mkind  = get_mkind(off_diag,map_operator_structure)
    pm = pauli_matrix[component]

    S = Operator()
    spin_range = range(len(spin_names))
    for n1, n2 in product(spin_range,spin_range):
        for on in orb_names:
            S += 0.5 * c_dag(*mkind(spin_names[n1],on)) * pm[n1,n2] * c(*mkind(spin_names[n2],on))
    return S

def S2_op(spin_names, orb_names, off_diag = None, map_operator_structure = None):
    Sz, Sp, Sm = map(lambda k: S_operator(k,spin_names,orb_names,off_diag,map_operator_structure), ('z','+','-'))
    return Sz*Sz + 0.5*(Sp*Sm + Sm*Sp)

def L_op(component, spin_names, orb_names, off_diag = None, map_operator_structure = None, basis='spherical', T=None):
    # FIXME
    assert component != 'y', "We cannot construct operators with complex coefficients at the moment. Sorry for that!"

    l = (len(orb_names)-1)/2
    L_melem_dict = {'z' : lambda m,mp: m if m==mp else 0,
                    '+' : lambda m,mp: np.sqrt(l*(l+1)-mp*(mp+1)) if m==mp+1 else 0,
                    '-' : lambda m,mp: np.sqrt(l*(l+1)-mp*(mp-1)) if m==mp-1 else 0,
                    'x' : lambda m,mp: 0.5*(L_melem_dict['+'](m,mp) + L_melem_dict['-'](m,mp)),
                    'y' : lambda m,mp: -0.5j*(L_melem_dict['+'](m,mp) - L_melem_dict['-'](m,mp))}
    L_matrix = np.zeros((2*l+1,2*l+1))
    orb_range = range(2*l+1)
    for o1, o2 in product(orb_range,orb_range): L_matrix[o1,o2] = L_melem_dict[component](o1-l,o2-l)

    # Transform from spherical basis if needed
    if basis == "cubic": T = spherical_to_cubic(l)
    if basis == "other" and T is None: raise ValueError("L_operator: provide T for other bases.")
    if T is not None: L_matrix = np.einsum("ij,jk,kl",np.conj(T),L_matrix,np.transpose(T))

    mkind = get_mkind(off_diag,map_operator_structure)
    L = Operator()
    for sn in spin_names:
        for o1, o2 in product(orb_range,orb_range):
            L += c_dag(*mkind(sn,orb_names[o1])) * L_matrix[o1,o2] * c(*mkind(sn,orb_names[o2]))
    return L

def L2_op(spin_names, orb_names, off_diag = None, map_operator_structure = None, basis='spherical', T=None):
    Lz, Lp, Lm = map(lambda k: L_operator(k,spin_names,orb_names,off_diag, map_operator_structure, basis, T),('z','+','-'))
    return Lz*Lz + 0.5*(Lp*Lm + Lm*Lp)

def LS_op(spin_names, orb_names, off_diag = None, map_operator_structure = None, basis='spherical', T=None):
    Sz, Sp, Sm = map(lambda k: S_operator(k,spin_names,orb_names,off_diag,map_operator_structure), ('z','+','-'))
    Lz, Lp, Lm = map(lambda k: L_operator(k,spin_names,orb_names,off_diag, map_operator_structure, basis, T),('z','+','-'))
    return Lz*Sz + 0.5*(Lp*Sm + Lm*Sp)