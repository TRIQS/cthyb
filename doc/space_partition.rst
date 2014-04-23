
Partitioning of the local Hilbert space
=======================================

The most time-consuming part of the Hybridization Expansion algorithm is a calculation
of dynamical traces, which eventually boils down to multiplication of matrices.
With the naïve approach it takes :math:`O(n^3)` operations to multiply two :math:`n\times n`
matrices (more sophisticated algorithms can reduce the complicity to approximately :math:`O(n^{2.4})`).

At the same time the dimension of the local Hilbert space grows exponentially with the number of
degrees of freedom, e.i. with the number of correlated orbitals and/or sites in the correlated cluster
under consideration.

It is then crucial to find a way to minimize the size of matrices being multiplied.
This goal can be achieved by partitioning the local Hilbert space of dimension :math:`N`
into a direct sum of :math:`K` subspaces with dimensions :math:`0 < N_k \le N`.
Let's say we have a basis set :math:`\{\psi_n\}` of the full Hilbert space,
which is to be partitioned into subsets :math:`\{\psi\}_1, \{\psi\}_2 \ldots \{\psi\}_K`.
These subsets are subject to two conditions:

#. The local Hamiltonian :math:`\mathcal{H}` maps each of the corresponding subspaces to itself:

    .. math::

        \forall k,\ \forall |\phi\rangle\in\mathrm{span}(\{\psi\}_k),\quad \mathcal{H}|\phi\rangle \in\mathrm{span}(\{\psi\}_k)

#. All creation and annihilation operators found within the dynamical traces map one subspace to one subspace or to zero.

    .. math::

        \forall \alpha,\ \forall k,\ \forall |\phi\rangle\in\mathrm{span}(\{\psi\}_k),\ \exists l,
        \quad c_\alpha |\phi\rangle \in\mathrm{span}(\{\psi\}_l)\\
        \forall \alpha,\ \forall k,\ \forall |\phi\rangle\in\mathrm{span}(\{\psi\}_k),\ \exists l,
        \quad c^\dagger_\alpha |\phi\rangle \in\mathrm{span}(\{\psi\}_l)

The problem can be slightly reformulated. One has to find a permutation of basis vectors, after which
(1) the local Hamiltonian is block-diagonal and (2) all :math:`c` and :math:`c^\dagger` operators
are block matrices with at most one non-zero block in each row and column. Such a permutation would group
basis states belonging to the same :math:`\{\psi\}_k` together.
        
Once these two requirements are met, the dynamical trace can be expressed as a sum of :math:`K` "strings".
Each string is a product of matrix blocks taken from :math:`\ c^\dagger,\ c` and :math:`\exp(-\Delta\tau\mathcal{H})`.
The complicity is reduced not only due to the smaller sizes of the blocks. In many cases
creation/annihilation operators map a subspace to zero (for example, annihilation of a particle in
a state with no particles). All strings including such mappings are immediately zero and can be
skipped altogether (the mappings are precomputed and stored at the very beginning of the simulation).

The present solver implements two partitioning strategies: partitioning with user-supplied quantum numbers
and the automatic partitioning. The latter strategy is more universal and is chosen by default.
One can override the choice by setting the parameter ``use_quantum_numbers`` to ``True``.

Quantum numbers
===============

This is the traditional approach to partitioning of the Hilbert space. User passes a list of Abelian
integrals of motion (operators) :math:`Q_1,\ldots,Q_L` as ``quantum_numbers`` parameter on construction
of the solver. Expectation values of these operators are calculated for each basis state :math:`\psi`, which
gives a combination of quantum numbers associated with the state:

    .. math::
        
        \langle\psi|Q_1|\psi\rangle, \ldots, \langle\psi|Q_L|\psi\rangle\ \Rightarrow q_1,\ldots, q_L
        
All states sharing the same combination of quantum numbers belong to the same subspace.

The first condition is obviously fulfilled by the obtained subspaces, because the Hamiltonian cannot
connect states with different quantum numbers by definition. If all operators :math:`Q_l` correspond
to Abelian symmetries of the system (like particle number and :math:`S_z`), the second condition is also
fulfilled.

This approach works well, but it requires some prior analysis of the local Hamiltonian from the user.
It may be difficult to discover an exhaustive set of integrals of motion, if the dimension of the local
Hilbert space is large and the interaction form is complicated.

Automatic partitioning
======================

The automatic partitioning algorithm employs no additional *a priori* information about the Hamiltonian
:math:`\mathcal{H}`. The only input data are the full set of basis states and the Hamiltonian itself.

The algorithm consists of two consequent steps (parts). On the first step it constructs the finest possible partition
which satisfies condition (1) alone. During the second step this partition is modified in a way to become
compatible also with condition (2).

**Part 1**

In the beginning the algorithm creates a data structure, which stores information about how the :math:`N` basis
vectors are partitioned into a number of subsets. Initially each basis vector resides alone in its own subset.

Then comes the main loop of the algorithm. The Hamiltonian is sequentially applied to each basis state
(initial state). Each application gives a linear combination of the basis states with only a few non-zero
coefficients, since H is normally sparse.
The algorithm iterates (inner loop) over all those basis vectors with non-zero coefficients (final states).

If the initial state and the final state reside in different subsets, these subsets are merged together.

Once the main loop is over, the partition of the basis is done. Two basis vectors are guaranteed to be in different subsets,
if they cannot be reached from each other by application of :math:`\mathcal{H}` any number of times.

All matrix elements of :math:`\mathcal{H}` are implicitly calculated along the line and stored for later use.

**Part 2**

In this part some subsets are additionally merged. The algorithm is applied in turn to all
:math:`c^\dagger_\alpha, c_\alpha` -pairs and for each pair it proceeds as follows:

1. Choose one of existing subsets, say :math:`\{\psi\}_m`.
2. Apply :math:`c_\alpha^\dagger` to all basis states in :math:`\{\psi\}_m`.
3. Iterate over all non-zero amplitudes of all resulting states and find all subsets, the corresponding basis states belong to.
   If no such amplitudes are present, skip to step 9.
4. Merge all found subsets to form a bigger subset :math:`\{\psi\}_{m'}`.
5. Apply :math:`c_\alpha` to all basis states in :math:`\{\psi\}_{m'}`.
6. Iterate over all non-zero amplitudes of all resulting states and find all subsets, the corresponding basis states belong to.
   If no such amplitudes are present, skip to step 9.
7. Merge all found subsets to form a bigger subset :math:`\{\psi\}_{m''}`.
8. Repeat steps 2-7 starting with :math:`\{\psi\}_{m''}`. Continue repeating those steps until no merges occur during the last iteration. 
9. Repeat steps 1-8 using another initial subset :math:`\{\psi\}_m`. If all subsets are exhausted, stop.

This procedure guarantees, that all creation and annihilation operators connect exactly one subspace to one subspace or to numerical zero.
The connections are stored (on step 8) and used afterwards to construct strings in the dynamical trace.

An important remark should be made about the data structure which maintains the partition of :math:`\{\psi\}`.
We use a *disjoint-set data structure* designed especially for fast ``find_set``, ``union_set`` and ``make_set``
operations [#ds]_, [#ds_boost]_. This structure is usually implemented as a forest with trees representing
disjoint subsets. Thanks to two special techniques called 'union by rank' and 'path compression', the amortized time
per operation is only :math:`O(\alpha(n))`. Here :math:`n` is the total number of ``find_set``, ``union_set`` and ``make_set``
operations, and :math:`\alpha(n)` is the extremely slowly-growing inverse Ackermann function. :math:`\alpha(n) < 5` for any
practical value of :math:`n`.

.. [#ds] R.E.Tarjan. *Data Structures and Network Algorithms.* Society for Industrial and Applied Mathematics, 1983.
.. [#ds_boost] Boost Disjoint Sets (http://www.boost.org/doc/libs/release/libs/disjoint_sets/disjoint_sets.html).
