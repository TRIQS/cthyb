#include "cthyb/solver_core.hpp"
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/gfs.hpp>
#include <triqs/test_tools/gfs.hpp>

using namespace cthyb;
using triqs::operators::c;
using triqs::operators::c_dag;
using triqs::operators::n;
using namespace triqs::gfs;
using triqs::hilbert_space::gf_struct_t;

TEST(CtHyb, Kanamori) {

  std::cout << "Welcome to the CTHYB solver\n";

  // Initialize mpi
  int rank = triqs::mpi::communicator().rank();

  // Parameters
  double beta      = 10.0;
  int num_orbitals = 2;
  double mu        = 1.0;
  double U         = 2.0;
  double J         = 0.2;
  double V         = 1.0;
  double epsilon   = 2.3;

  auto N     = [](std::string sn, int an) { return n(sn + '-' + std::to_string(an), 0); };
  auto C     = [](std::string sn, int an) { return c(sn + '-' + std::to_string(an), 0); };
  auto C_dag = [](std::string sn, int an) { return c_dag(sn + '-' + std::to_string(an), 0); };

  // Hamiltonian
  many_body_op_t H;
  for (int o = 0; o < num_orbitals; ++o) { H += U * N("up", o) * N("down", o); }
  for (int o1 = 0; o1 < num_orbitals; ++o1)
    for (int o2 = 0; o2 < num_orbitals; ++o2) {
      if (o1 == o2) continue;
      H += (U - 2 * J) * N("up", o1) * N("down", o2);
    }
  for (int o1 = 0; o1 < num_orbitals; ++o1)
    for (int o2 = 0; o2 < num_orbitals; ++o2) {
      if (o2 >= o1) continue;
      H += (U - 3 * J) * N("up", o1) * N("up", o2);
      H += (U - 3 * J) * N("down", o1) * N("down", o2);
    }

  for (int o1 = 0; o1 < num_orbitals; ++o1)
    for (int o2 = 0; o2 < num_orbitals; ++o2) {
      if (o1 == o2) continue;
      H += -J * C_dag("up", o1) * C_dag("down", o1) * C("up", o2) * C("down", o2);
      H += -J * C_dag("up", o1) * C_dag("down", o2) * C("up", o2) * C("down", o1);
    }

#ifdef QN
  // quantum numbers
  std::vector<many_body_op_t> qn;
  qn.resize(2);
  for (int o = 0; o < num_orbitals; ++o) {
    qn[0] += N("up", o);
    qn[1] += N("down", o);
  }
#endif

  // gf structure
  gf_struct_t gf_struct;
  for (int o = 0; o < num_orbitals; ++o) gf_struct.push_back({"down-" + std::to_string(o),{0}});
  for (int o = 0; o < num_orbitals; ++o) gf_struct.push_back({"up-" + std::to_string(o),{0}});

  // Construct CTQMC solver
  solver_core solver({beta, gf_struct, 1025, 2500});

  // Set G0
  triqs::clef::placeholder<0> om_;
  auto g0_iw = gf<imfreq>{{beta, Fermion}, {1, 1}};
  g0_iw(om_) << om_ + mu - (V / (om_ - epsilon) + V / (om_ + epsilon));
  for (int o = 0; o < 2 * num_orbitals; ++o) { solver.G0_iw()[o] = triqs::gfs::inverse(g0_iw); }

  // Solve parameters
  auto n_cycles     = 5000;
  auto p            = solve_parameters_t(H, n_cycles);
  p.max_time        = -1;
  p.random_name     = "";
  p.random_seed     = 123 * rank + 567;
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
  std::string filename = "kanamori";
#ifdef QN
  filename += "_qn";
#endif

  auto & G_tau = *solver.G_tau;

  if (rank == 0) {
    triqs::h5::file G_file(filename + ".out.h5", 'w');
    for (int o = 0; o < num_orbitals; ++o) {
      h5_write(G_file, "G_up-" + std::to_string(o), G_tau[o]);
      h5_write(G_file, "G_down-" + std::to_string(o), G_tau[num_orbitals + o]);
    }
  }

  gf<imtime> g;
  if (rank == 0) {
    triqs::h5::file G_file(filename + ".ref.h5", 'r');
    for (int o = 0; o < num_orbitals; ++o) {
      h5_read(G_file, "G_up-" + std::to_string(o), g);
      EXPECT_GF_NEAR(g, G_tau[o]);
      h5_read(G_file, "G_down-" + std::to_string(o), g);
      EXPECT_GF_NEAR(g, G_tau[num_orbitals + o]);
    }
  }
}
MAKE_MAIN;
