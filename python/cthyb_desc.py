from wrap_generator import *

# The cthyb module
module = module_(full_name = "pytriqs.applications.impurity_solvers.cthyb", doc = "The cthyb matrix solver")
module.use_module('gf')
module.use_module('parameters')
module.use_module('operators2')
module.add_include("c++/solver_core.hpp")
module.add_using("namespace cthyb")

# The Solver class
cthyb = class_(
        py_type = "SolverCore",
        c_type = "solver_core"
        )

module.add_class(cthyb)

cthyb.add_constructor(signature = "(double beta, std::map<std::string, std::vector<int>> gf_struct, int n_iw=5000, int n_tau=10001, int n_l=50)", doc = """ """)

cthyb.add_method(name = "solve",
             signature = "void(triqs::utility::many_body_operator<double> h_loc, params::parameters params, std::vector<triqs::utility::many_body_operator<double>> quantum_numbers = std::vector<triqs::utility::many_body_operator<double>>{}, bool use_quantum_numbers = false)",
             doc = """ """)

cthyb.add_method(name = "solve_parameters",
                 signature = "parameters()",
                 is_static = True,
                 doc = """Get the form of solve parameters""")

cthyb.add_property(name = "G0_iw",
                   getter = cfunction("block_gf_view<imfreq> G0_iw()"),
                   setter = cfunction("void set_G0_iw(block_gf_view<imfreq> G0)"),
                   doc = "G0(iw) in imaginary frequencies")

cthyb.add_property(name = "Delta_tau",
                   getter = cfunction("block_gf_view<imtime> Delta_tau()"),
                   doc = "Delta(tau) in imaginary time")

cthyb.add_property(name = "G_tau",
                   getter = cfunction("block_gf_view<imtime> G_tau()"),
                   doc = "G(tau) in imaginary time")

cthyb.add_property(name = "G_l",
                   getter = cfunction("block_gf_view<legendre> G_l()"),
                   doc = "G_l in Legendre polynomials representation")

cthyb.add_property(name = "atomic_gf",
                   getter = cfunction("block_gf_view<imtime> atomic_gf()"),
                   doc = "Atomic G(tau) in imaginary time")

module.add_function("gf_view<imtime> change_mesh(gf_view<imtime> old_gf, int new_n_tau)")
module.add_function("block_gf_view<imtime> change_mesh(block_gf_view<imtime> old_gf, int new_n_tau)")

# generate the module code
module.generate_code()
