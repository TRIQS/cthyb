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

TEST(CtHyb, Legendre) {

  std::cout << "Welcome to the CTHYB solver\n";

  // Initialize mpi
  int rank = triqs::mpi::communicator().rank();

  // Parameters
  double beta     = 10.0;
  double U        = 2.0;
  double mu       = 1.0;
  double h        = 0.0;
  double V        = 1.0;
  double epsilon1 = 2.1;
  double epsilon2 = -2.4;

  // define operators and QN
  auto H = U * n("up", 0) * n("down", 0) + h * n("up", 0) - h * n("down", 0);
  // quantum numbers
  std::vector<many_body_op_t> qn{n("up", 0), n("down", 0)};
  // gf structure
  std::map<std::string, indices_type> gf_struct{{"up", {0}}, {"down", {0}}};

  // Construct CTQMC solver
  solver_core solver(beta, gf_struct, 1025, 2500, 50);

  // Set G0
  triqs::clef::placeholder<0> om_;
  auto g0_iw = gf<imfreq>{{beta, Fermion}, {1, 1}};
  g0_iw(om_) << om_ + mu - (V * V / (om_ - epsilon1) + V * V / (om_ - epsilon2));
  for (int bl = 0; bl < 2; ++bl) solver.G0_iw()[bl] = triqs::gfs::inverse(g0_iw);

  // Solve parameters
  int n_cycles       = 5000;
  auto p             = solve_parameters_t(H, n_cycles);
  p.random_name      = "";
  p.random_seed      = 123 * rank + 567;
  p.max_time         = -1;
  p.length_cycle     = 50;
  p.n_warmup_cycles  = 50;
  p.measure_g_tau    = false;
  p.measure_g_l      = true;
  p.quantum_numbers  = qn;
  p.partition_method = "quantum_numbers";
  p.move_double      = false;

  // Solve!
  solver.solve(p);

  // Save the results
  if (rank == 0) {
    triqs::h5::file G_file("legendre.out.h5", 'w');
    h5_write(G_file, "G_up", solver.G_l()[0]);
    h5_write(G_file, "G_down", solver.G_l()[1]);
  }

  gf<legendre> g;
  if (rank == 0) {
    triqs::h5::file G_file("legendre.ref.h5", 'r');
    h5_read(G_file, "G_up", g);
    EXPECT_GF_NEAR(g, solver.G_l()[0]);
    h5_read(G_file, "G_down", g);
    EXPECT_GF_NEAR(g, solver.G_l()[1]);
  }
}
MAKE_MAIN;
