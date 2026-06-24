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
   guide/CRM_Dyson_solver


Tutorials
---------

.. toctree::
   :maxdepth: 1

   guide/aim
   guide/slater_five_band
   guide/cthyb_convergence_tests
   guide/high_freq_moments

Reference manual
----------------

.. autosummary::
   :toctree: _ref
   :template: autosummary_module_template.rst
   :recursive:

   triqs_cthyb.configuration
   triqs_cthyb.multiplet_tools
   triqs_cthyb.solver
   triqs_cthyb.solver_core
   triqs_cthyb.tail_fit
   triqs_cthyb.util

Links to all relevant solver parameters:

* :doc:`Construction parameters <_ref/triqs_cthyb.solver_core.ConstrParametersT>`
* :doc:`Solve parameters <_ref/triqs_cthyb.solver_core.SolveParametersT>`

The C++ reference manual can be found `here <./doxygen/index.html>`_.

FAQs
----

.. toctree::
   :maxdepth: 2

   faqs
