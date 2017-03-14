#include "solver_core.hpp"
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/gfs.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace cthyb;
using triqs::operators::c;
using triqs::operators::c_dag;
using triqs::operators::n;
using namespace triqs::gfs;
using indices_type = triqs::operators::indices_t;

TEST(CtHyb, AtomicGf) {

  std::cout << "Welcome to the CTHYB solver\n";

  // Initialize mpi
  int rank = triqs::mpi::communicator().rank();

  // Parameters
  double beta = 10.0;
  double U0   = 2.0;
  double U1   = 2.0;
  double ed0  = -1.1;
  double ed1  = -0.9;
  double V    = 1.5;

  // define operators
  auto H = U0 * n("up", 0) * n("dn", 0) + U1 * n("up", 1) * n("dn", 1);
  H += ed0 * (n("up", 0) + n("dn", 0)) + ed1 * (n("up", 1) + n("dn", 1));
  H += V * (c_dag("up", 0) * c("up", 1) + c_dag("up", 1) * c("up", 0) + c_dag("dn", 0) * c("dn", 1) + c_dag("dn", 1) * c("dn", 0));

  std::map<std::string, indices_type> gf_struct{{"up", {0, 1}}, {"dn", {0, 1}}};

  // Construct CTQMC solver
  solver_core solver(beta, gf_struct, 1025, 2051);

  // Solve parameters
  auto p            = solve_parameters_t(H, 0);
  p.length_cycle    = 1;
  p.n_warmup_cycles = 0;

  triqs::clef::placeholder<0> om_;
  auto g0_iw = gf<imfreq>{{beta, Fermion}, {2, 2}};
  g0_iw(om_) << om_;
  solver.G0_iw()[0] = triqs::gfs::inverse(g0_iw);
  solver.G0_iw()[1] = triqs::gfs::inverse(g0_iw);

  // Solve!
  solver.solve(p);

  // Save the results
  if (rank == 0) {
    triqs::h5::file G_file("atomic_gf.out.h5", 'w');
    h5_write(G_file, "G_tau", solver.atomic_gf());
  }

  block_gf<imtime> g;
  if (rank == 0) {
    triqs::h5::file G_file("atomic_gf.ref.h5", 'r');
    h5_read(G_file, "G_tau", g);
    EXPECT_BLOCK_GF_NEAR(g, solver.atomic_gf());
  }
}
MAKE_MAIN;
