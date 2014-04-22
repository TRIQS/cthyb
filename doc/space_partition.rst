
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

        \forall k,\ \forall \phi\in\mathrm{span}(\{\psi\}_k),\quad \mathcal{H}\phi \in\mathrm{span}(\{\psi\}_k)

#. All creation and annihilation operators found within the dynamical traces map one subspace to one subspace or to zero.

    .. math::

        \forall \alpha,\ \forall k,\ \forall \phi\in\mathrm{span}(\{\psi\}_k),\ \exists l,
        \quad c_\alpha \phi \in\mathrm{span}(\{\psi\}_l)\\
        \forall \alpha,\ \forall k,\ \forall \phi\in\mathrm{span}(\{\psi\}_k),\ \exists l,
        \quad c^\dagger_\alpha \phi \in\mathrm{span}(\{\psi\}_l)

The problem can be slightly reformulated. One has to find a permutation of basis vectors, after which
1) the local Hamiltonian is block-diagonal and 2) all :math:`c` and :math:`c^\dagger` operators
are block matrices with at most one non-zero block in each row and column. Such a permutation would group
basis states belonging to the same :math:`\{\psi\}_k` together.
        
Once these two requirements are met, the dynamical trace can be expressed as a sum of :math:`K` "strings".
Each string is a product of matrix blocks taken from :math:`\ c^\dagger,\ c` and :math:`\exp(-\Delta\tau\mathcal{H})`.
The complicity is reduced not only due to the smaller sizes of the blocks. In many cases
creation/annihilation operators map a subspace to zero (for example, annihilation of a particle in
a state with no particles). All strings including such mappings are immediately zero and can be
skipped altogether (the mappings are precomputed and stored at the very beginning of the simulation).

This solver implements two partitioning strategies: partitioning with user-supplied quantum numbers
and the automatic partitioning. The latter strategy is more universal and is chosen by default.
One can override the choice by setting the parameter ``use_quantum_numbers`` to ``True``.

Quantum numbers
===============

This is the traditional approach to partitioning of the Hilbert space. User passes a list of abelian
integrals of motion (operators) :math:`Q_1,\ldots,Q_L` as ``quantum_numbers`` parameter on construction
of the solver. Expectation values of these operators are calculated for each basis state :math:`\psi`, which
gives a combination of quantum numbers associated with the state:

    .. math::
        
        \langle\psi|Q_1|\psi\rangle, \ldots, \langle\psi|Q_L|\psi\rangle\ \Rightarrow q_1,\ldots, q_L
        
All states sharing the same combination of quantum numbers belong to the same subspace.

The first condition is obviously fulfilled by the obtained subspaces, because the Hamiltonian cannot
connect states with different quantum numbers by definition. If all operators :math:`Q_l` correspond
to abelian symmetries of the system (like particle number and :math:`S_z`), the second condition is also
fulfilled.

This approach works well, but it requires some prior analysis of the local Hamiltonian from the user.
It may be difficult to discover an exhaustive set of integrals of motion, if the dimension of the local
Hilbert space is large and the interaction form is complicated.