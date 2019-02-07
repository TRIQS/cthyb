.. _moves:

QMC moves
=========

In ``cthyb``, a *configuration* (state of the Markov chain) is a :math:`\tau`-ordered sequence
of an even number of creation and annihilation operators :math:`c^\dagger_{Ai}(\tau), c_{Bj}(\tau)` with block
indices :math:`A, B` and inner indices :math:`i, j`.

Here is a list of all moves (updates) used to organize random walk in the space of possible configurations.

Insert one pair of operators
****************************

Randomly choose a block index :math:`A`, two inner indices :math:`i, j` from the block :math:`A`, and
two imaginary time points :math:`\tau, \tau'`. Try to insert a pair :math:`c^\dagger_{Ai}(\tau), c_{Aj}(\tau')`
into the configuration.

Probability of choosing a particular block index :math:`A` can be adjusted through the parameter ``proposal_prob``.

This move is always enabled.

Remove one pair of operators
****************************

Randomly choose a block index :math:`A`, and two operators with this block index :math:`c^\dagger_{Ai}(\tau), c_{Aj}(\tau')`
from the configuration. Try to remove the chosen operators.

Probability of choosing a particular block index :math:`A` can be adjusted through the parameter ``proposal_prob``.

This move is always enabled.

Insert two pairs of operators
*****************************

Similarly to the one pair insertion, randomly choose two combinations :math:`(A;i_1,j_1;\tau_1,\tau'_1)` and
:math:`(B;i_2,j_2;\tau_2,\tau'_2)`. Try to insert two pairs of operators
:math:`c^\dagger_{Ai_1}(\tau_1), c_{Aj_1}(\tau'_1), c^\dagger_{Bi_2}(\tau_2), c_{Aj_2}(\tau'_2),` into the configuration.

Probability of choosing a particular block index can be adjusted through the parameter ``proposal_prob``.
:math:`A` and :math:`B` are chosen independently.

This move is enabled by default and can be disabled by setting ``move_double`` to ``False``.

Remove two pairs of operators
*****************************

Randomly choose two block indices :math:`A` and :math:`B`, two operators with the block index
:math:`A`, :math:`c^\dagger_{Ai_1}(\tau_1), c_{Aj_1}(\tau'_1)` and two operators with the block index
:math:`B`, :math:`c^\dagger_{Bi_2}(\tau_2), c_{Bj_2}(\tau'_2)`. Try to remove the chosen operators.

Probability of choosing a particular block index can be adjusted through the parameter ``proposal_prob``.
:math:`A` and :math:`B` are chosen independently.

.. note::

    The two-pairs insertion/removal moves are **mandatory** for multiorbital models with interaction
    Hamiltonians beyond the density-density approximation! It was shown that the single-pair moves alone
    cannot reach some specific configurations in these cases. Those configurations are physically relevant and
    contribute to the observables.

Shift one operator
******************

Randomly choose one operator from the configuration. Choose a random inner index from the same block and
a random time point. Try to move the chosen operator to the new time position and replace its inner index.

This move helps to reduce statistical noise. It is enabled by default and can be disabled with ``move_shift = False``.

Global move - change of operator indices
****************************************

Randomly choose one of the user-provided index substitution maps,
`(old block index, old inner index) -> (new block index, new inner index)`.
Find all operators in the configuration, whose indices are to be changed according to the
chosen substitution map. Let the total number of such operators be :math:`N_{max}`.
A uniformly distributed integer :math:`N\in[1;N_{max}]` is then generated, and
:math:`N` randomly picked operators (out of :math:`N_{max}`) actually change their indices.

This move is disabled by default, because it is much slower than the rest of the moves.
It can alleviate some ergodicity problems occuring at lower temperatures (confinement in a local minimum).
*However, one may not assume that all problems of this kind are automatically solved by the move.*

Proposal probability of the global move is set by ``move_global_prob`` parameter.
The substitution maps must be passed as a dictionary to ``move_global`` parameter.
The dictionary has the following format:

``dict(map_name : dict((A_1,i_1) : (B_1,j_1), (A_2,i_2) : (B_2,j_2), ...), ...)``

An empty dictionary (default) disables the move completely.
