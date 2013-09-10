
Implementation notes
====================

How to read the code
--------------------

The class containing the solver is called ``ctqmc_seg`` and is defined in

* :file:`ctqmc_seg.hpp`

The solver is built on a parameter dictionary. In order to run, it needs moves
and measures, which are defined in the files starting with :file:`move_` and
:file:`measure_`, respectively:

* :file:`move_insert_segment.hpp`
* :file:`move_remove_segment.hpp`
* :file:`measure_nn.hpp`

Note that so far I haven't added other measures.  The moves and measures act
on, or use, a set of data characterizing the Monte-Carlo process: this set of
data is defined in

* :file:`qmc_data.hpp`

A ``qmc_data`` object contains all the data that the moves and measures need, including cached data:

* input parameters (like the inverse temperature, the interaction matrix U, the chemical potential mu...)
* the configuration (=the segments on their respective lines)
* cached information (like the determinants, the overlap matrix and the sign)
 
In principle, the last item could be seen as part of the configuration. A
configuration, however, is entirely *defined* by the locations of its segments'
start and end points and full/empty line configuration. From this mere
information, one could recompute the determinants, the overlap matrix, the
sign, as well as other observables that one may want to cache to gain speed.


Configuration
-------------

The configuration is defined in

* :file:`configuration.hpp`

It contains helper methods to manipulate configurations of segments: 

* ``insert_segment``  --> segment insertion
* ``remove_segment``  --> segment removal
* ``find_segment``    --> find a segment given its index
* ``length_overlaps`` --> compute overlaps/length of a segment in a configuration
* ``maximal_length``  --> compute maximal length of a segment given an time position
* ``length``          --> segment length 

All these methods are used in moves to manipulate the configuration.  The
configuration also contains the internal storage of the segments, which is
defined in

* :file:`operator_maps.hpp`

This last class contains the low-level methods to add and remove operators in a
configuration. It is never accessed by other classes than the configuration
class. It mainly consists in two maps of operators:

* ``fullopmap`` --> contains all operators from all flavors
* ``opmaps[flavor]`` --> contains operators of the flavor *flavor*, pointing to the operators of the fullopmap.

Iterators on these two maps have been wrapped to simplify access to their members:

* ``flavored_const_iterator`` is used to iterate on a given flavor's line
* ``const_iterator`` is used to iterate among iterators, irrespective of their flavors.

Both have a method ``cyclic_right()`` and ``cyclic_left()`` which allow to iterator
through operators **with periodicity**. Dereferencing these iterators yield
information about the operator:

* ``it->tau``               --> imaginary time position
* ``it->dagger``            --> creation/annihilation operator
* ``it->flavor``            --> flavor of the operator
* ``it->right_occupations`` --> occupations of all lines right of the operator


