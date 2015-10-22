.. _slater:

Another example: a multiorbital impurity model
==============================================

We present another more complex example here. This time, we solve the five-band
impurity with a fully-rotationally invariant interaction:

.. math::

  \mathcal{H}_\mathrm{int} = \frac{1}{2} \sum_{ijkl,\sigma \sigma'} U_{ijkl} a_{i \sigma}^\dagger a_{j \sigma'}^\dagger a_{l \sigma'} a_{k \sigma}.

Here is the python :download:`script <slater_five_band.py>`:

.. literalinclude:: slater_five_band.py

This script generates an HDF5 archive file called :file:`slater_five_band.h5`.
This file contains the Green's function in imaginary time and in imaginary
frequencies found by the solver. Let us plot the Green's function:

.. plot:: guide/slater_five_band_plot.py
   :include-source:
   :scale: 70

We have rebinned the data to 1000 points to reduce the noise. The calculation
takes about an hour and data was accumulated on 32 cores.
