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

TEST(CtHyb, Anderson) {

  std::cout << "Welcome to the CTHYB solver\n";

  // Initialize mpi
  int rank = triqs::mpi::communicator().rank();

  // Parameters
  double beta    = 10.0;
  double U       = 2.0;
  double mu      = 1.0;
  double h       = 0.0;
  double V       = 1.0;
  double epsilon = 2.3;

  // GF structure
  enum spin { up, down };
#ifdef BLOCK
  std::map<std::string, indices_type> gf_struct{{"up", {0}}, {"down", {0}}};
  auto n_up   = n("up", 0);
  auto n_down = n("down", 0);
#else
  std::map<std::string, indices_type> gf_struct{{"tot", {0, 1}}};
  auto n_up   = n("tot", 0);
  auto n_down = n("tot", 1);
#endif

  // define operators
  auto H = U * n_up * n_down + h * n_up - h * n_down;

#ifdef QN
  // Quantum numbers
  std::vector<many_body_op_t> qn{n_up, n_down};
#endif

  // Construct CTQMC solver
  solver_core solver(beta, gf_struct, 1025, 2500);

  // Set G0
  triqs::clef::placeholder<0> om_;
#ifdef BLOCK
  auto g0_iw = gf<imfreq>{{beta, Fermion}, {1, 1}};
  g0_iw(om_) << om_ + mu - (V * V / (om_ - epsilon) + V * V / (om_ + epsilon));
  for (int bl = 0; bl < 2; ++bl) solver.G0_iw()[bl] = triqs::gfs::inverse(g0_iw);
#else
  auto g0_iw = gf<imfreq>{{beta, Fermion}, {2, 2}};
  g0_iw(om_) << om_ + mu - (V * V / (om_ - epsilon) + V * V / (om_ + epsilon));
  solver.G0_iw()[0] = triqs::gfs::inverse(g0_iw);
#endif

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
  std::string filename = "anderson";
#ifdef BLOCK
  filename += "_block";
#endif
#ifdef QN
  filename += "_qn";
#endif

  if (rank == 0) {
    triqs::h5::file G_file(filename + ".out.h5", 'w');
#ifdef BLOCK
    h5_write(G_file, "G_up", solver.G_tau->operator[](0));
    h5_write(G_file, "G_down", solver.G_tau->operator[](1));
#else
    h5_write(G_file, "G_tot", solver.G_tau->operator[](0));
#endif
  }

  gf<imtime> g;
  if (rank == 0) {
    triqs::h5::file G_file(filename + ".ref.h5", 'r');
#ifdef BLOCK
    h5_read(G_file, "G_up", g);
    EXPECT_GF_NEAR(g, solver.G_tau->operator[](0));
    h5_read(G_file, "G_down", g);
    EXPECT_GF_NEAR(g, solver.G_tau->operator[](1));
#else
    h5_read(G_file, "G_tot", g);
    EXPECT_GF_NEAR(g, solver.G_tau->operator[](0));
#endif
  }
}
MAKE_MAIN;
