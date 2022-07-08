
.. _faqs:

Frequently-asked questions
==========================

Q: Why is my code running so slowly?
------------------------------------

A: If you see massive performance problems (slow down by huge factors like
100), then you might have an issue with OpenMP when using MKL for instance.
Setting ``OMP_NUM_THREADS=1`` solves this problem.


Q: How do I save the triqs and cthyb hash and script for debugging purposes?
----------------------------------------------------------------------------

A: Simply add this to your script::

    from triqs_cthyb import version
    if mpi.is_master_node():
      with HDFArchive(filename+".h5",'a') as Results:
        if "log" not in Results: Results.create_group("log")
        log = Results["log"]
        log["version"] = version.version
        log["triqs_hash"] = version.triqs_hash
        log["cthyb_hash"] = version.cthyb_hash
        log["script"] = open(sys.argv[0]).read() # read myself !

Q: Why does my data look so noisy?
----------------------------------

A: If you are running a parallel calculation, ensure that you are using a
different random seed on each core, i.e., that it is a function of the MPI
rank::

    param['random_seed'] = 34788 + 928374 * mpi.rank   # Default random seed

Q: How do I use the segment picture?
------------------------------------

A: This cthyb code is based on the matrix formulation and does not include
optimisations for the segment picture (applicable only in cases with
density-density interactions only).
