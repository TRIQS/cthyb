#include <triqs_cthyb/solver_core.hpp>

#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace triqs_cthyb;
using triqs::operators::c;
using triqs::operators::c_dag;
using triqs::operators::n;
using namespace triqs::gfs;
using namespace triqs::mesh;
using triqs::hilbert_space::gf_struct_t;

TEST(CtHyb, Spinless) {

  std::cout << "Welcome to the CTHYB solver\n";

  // Initialize mpi
  int rank = mpi::communicator().rank();

  // Parameters
  double beta    = 10.0;
  double U       = 2.0;
  double mu      = 1.0;
  double epsilon = 2.3;
  double t       = 0.1;

  // define operators
  auto H = U * n("tot", 0) * n("tot", 1) - t * c_dag("tot", 0) * c("tot", 1) - t * c_dag("tot", 1) * c("tot", 0);
  gf_struct_t gf_struct{{"tot", 2}};

#ifdef QN
  // quantum numbers
  std::vector<many_body_op_t> qn{n("tot", 0) + n("tot", 1)};
#endif

  // Construct CTQMC solver
  solver_core solver({beta, gf_struct, 1025, 2500});

  // Set hybridization function
  nda::clef::placeholder<0> om_;
  auto delta_iw = gf<imfreq>{{beta, Fermion}, {2, 2}};

  using nda::range;
  auto d00 = slice_target(delta_iw, range(0, 1), range(0, 1));
  auto d11 = slice_target(delta_iw, range(1, 2), range(1, 2));
  auto d01 = slice_target(delta_iw, range(0, 1), range(1, 2));
  auto d10 = slice_target(delta_iw, range(1, 2), range(0, 1));
  d00(om_) << (om_ - epsilon) * (1.0 / (om_ - epsilon - t)) * (1.0 / (om_ - epsilon + t))
        + (om_ + epsilon) * (1.0 / (om_ + epsilon - t)) * (1.0 / (om_ + epsilon + t));
  d11(om_) << d00(om_);
  d01(om_) << -t * (1.0 / (om_ - epsilon - t)) * (1.0 / (om_ - epsilon + t)) - t * (1.0 / (om_ + epsilon - t)) * (1.0 / (om_ + epsilon + t));
  d10(om_) << d01(om_);

  // Set G0
  auto g0_iw = gf<imfreq>{{beta, Fermion}, {2, 2}};
  g0_iw(om_) << om_ + mu - delta_iw(om_);
  solver.G0_iw()[0] = triqs::gfs::inverse(g0_iw);

  // Solve parameters
  int n_cycles      = 5000;
  auto p            = solve_parameters_t(H, n_cycles);
  p.random_name     = "";
  p.random_seed     = 123 * rank + 567;
  p.max_time        = -1;
  p.length_cycle    = 50;
  p.n_warmup_cycles = 50;
  p.move_double     = false;
#ifdef QN
  p.quantum_numbers  = qn;
  p.partition_method = "quantum_numbers";
#endif

  // Solve!
  solver.solve(p);

  // Save the results
  std::string filename = "spinless";
#ifdef QN
  filename += "_qn";
#endif

  auto & G_tau = *solver.G_tau;

  if (rank == 0) {
    h5::file G_file(filename + ".out.h5", 'w');
    h5_write(G_file, "G_tau", G_tau[0]);
  }

  gf<imtime> g;
  if (rank == 0) {
    h5::file G_file(filename + ".ref.h5", 'r');
    h5_read(G_file, "G_tau", g);
    EXPECT_GF_NEAR(g, G_tau[0]);
  }
}
MAKE_MAIN;
