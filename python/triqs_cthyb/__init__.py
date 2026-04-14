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
DOC

"""
# triqs::stat::histogram so that Solver.perturbation_order(_total) and performance_analysis can be used
from triqs.stat.histograms import Histogram
from . import variants
from .solver import Solver
from .util import estimate_nfft_buf_size
from .solve_generic import solve_generic, TailFitParams, LegendreParams, CRMParams

# Importing an extension module fixes the solver variant of the process, so the names it
# provides are resolved on demand and kept out of __all__ (see variants.py).
_lazy_attrs = {'SolverCore': 'solver_core', 'ConstrParametersT': 'solver_core', 'SolveParametersT': 'solver_core',
               'Configuration': 'configuration', 'OpDesc': 'configuration'}

__all__ = ['Solver', 'estimate_nfft_buf_size',
           'solve_generic', 'TailFitParams', 'LegendreParams', 'CRMParams']


def __getattr__(name):
    module_name = _lazy_attrs.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    return getattr(variants.load(module_name), name)


def __dir__():
    return sorted(__all__ + list(_lazy_attrs))
