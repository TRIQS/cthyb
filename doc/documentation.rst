.. _documentation:

Documentation
*************

Basic notions
-------------

.. toctree::
   :maxdepth: 1

   basicnotions/cthyb
   basicnotions/moves
   basicnotions/measurements

User guide
----------

.. toctree::
   :maxdepth: 1
   
   guide/settingparameters
   guide/dmft
   guide/random
   guide/static_observables_notebook
   guide/multiplet_analysis_notebook
   guide/dynamic_susceptibility_notebook
   guide/perturbation_order_notebook
   

Tutorials
---------

.. toctree::
   :maxdepth: 1

   guide/aim
   guide/slater_five_band
   guide/cthyb_convergence_tests

Reference manual
----------------

.. autosummary::
   :toctree: _ref
   :template: autosummary_module_template.rst
   :recursive:

   triqs_cthyb.multiplet_tools
   triqs_cthyb.solver
   triqs_cthyb.tail_fit
   triqs_cthyb.util

Link to all relevant solver parameters:

.. toctree::
   :maxdepth: 1

   _ref/triqs_cthyb.solver.Solver.solve_parameters.rst
   _ref/triqs_cthyb.solver.Solver.constr_parameters.rst

FAQs
----

.. toctree::
   :maxdepth: 2

   faqs
