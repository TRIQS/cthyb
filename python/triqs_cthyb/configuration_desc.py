from cpp2py.wrap_generator import *

module = module_(full_name = "configuration", doc = r"The TRIQS cthyb configuration", app_name = "triqs_cthyb")

module.add_include("triqs_cthyb/solver_core.hpp")

# The class configuration
c = class_(
        py_type = "Configuration",  # name of the python class
        c_type = "triqs_cthyb::configuration",   # name of the C++ class
        doc = r"""Core class of the cthyb configuration""",   # doc of the C++ class
        comparisons = "==",
        is_printable = True,
        hdf5 = True
)

c.add_property(name = "beta",
               getter = cfunction("double beta()"),
               doc = r"""Value of beta""")

c.add_property(name = "size",
               getter = cfunction("int size()"),
               doc = r"""Size of the configuration""")

module.add_class(c)

module.generate_code()
