################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2025 by The Simons Foundation
#     author: N. Wentzell
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
Selection of the compiled solver variant.

The variants are separate extension modules that wrap the same C++ types. c2py
registers wrapped types in a table shared by every module of the process, keyed by
the C++ type, so the second variant to be imported looks up the Python types of the
first one and rejects its own objects. Only one variant can therefore be used per
process, and the extension modules are imported on demand so that the choice is not
made behind the back of ``Solver(hybridisation_is_complex=...)``.
"""

import importlib

# (hybridisation_is_complex, local_hamiltonian_is_complex) -> (subpackage, cmake option)
_variants = {
    (False, False): (None, None),
    (True, False): ("_complex_hyb", "-DHybridisation_is_complex=ON"),
    (True, True): ("_complex_all", "-DLocal_hamiltonian_is_complex=ON"),
}

# Extension modules that must be imported alongside another one, because it exchanges
# objects of the types they wrap.
_prerequisites = {"solver_core": ("configuration",)}

_selected = None


def _name(key):
    return _variants[key][0] or "real"


def select(hybridisation_is_complex=False, local_hamiltonian_is_complex=False):
    """Return the requested variant, fixing the variant of this process on first use."""

    key = (bool(hybridisation_is_complex), bool(local_hamiltonian_is_complex))
    if key not in _variants:
        raise ValueError("local_hamiltonian_is_complex requires hybridisation_is_complex=True")

    global _selected
    if _selected is None:
        _selected = key
    elif _selected != key:
        raise RuntimeError(
            f"The {_name(_selected)} solver variant is already loaded in this process, so the "
            f"{_name(key)} variant cannot be used as well. The variants wrap the same C++ types "
            "and only one of them can be imported per process. Use the other variant in a "
            "separate process, and import Solver rather than SolverCore, which loads the real "
            "variant.")
    return _selected


def load(module_name, key=None):
    """Import an extension module of the selected variant, and the ones it depends on."""

    if key is None:
        key = select()

    for prerequisite in _prerequisites.get(module_name, ()):
        load(prerequisite, key)

    subpackage, option = _variants[key]
    if subpackage is None:
        return importlib.import_module(f".{module_name}", __package__)

    try:
        return importlib.import_module(f".{subpackage}.{module_name}", __package__)
    except ImportError:
        raise ImportError(f"The {_name(key)} solver variant is not available. Rebuild with {option}") from None
