
.. _faqs:

Frequently-asked questions
==========================

Q: How do I save the triqs and cthyb hash and script for debugging purposes?
----------------------------------------------------------------------------

A: Simply add this to your script::

    from pytriqs.applications.impurity_solvers.cthyb import version
    if mpi.is_master_node():
      with HDFArchive(filename+".h5",'a') as Results:
        if "log" not in Results: Results.create_group("log")
        log = Results["log"]
        log["version"] = version.release
        log["triqs_hash"] = version.triqs_hash
        log["cthyb_hash"] = version.cthyb_hash
        log["script"] = open(sys.argv[0]).read() # read myself !
