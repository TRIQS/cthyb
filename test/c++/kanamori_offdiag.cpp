#include <triqs_cthyb/solver_core.hpp>

#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/test_tools/gfs.hpp>
#include <triqs/utility/itertools.hpp>

using namespace triqs_cthyb;
using namespace triqs::gfs;
using namespace triqs::mesh;
using triqs::operators::c;
using triqs::operators::c_dag;
using triqs::operators::n;
using triqs::arrays::matrix;
using triqs::hilbert_space::gf_struct_t;

TEST(CtHyb, KanamoriOffDiag) {

  std::cout << "Welcome to the CTHYB solver\n";

  // Initialize mpi
  int rank = mpi::communicator().rank();

  // Parameters
  double beta            = 10.0;
  const int num_orbitals = 2;
  double mu              = 1.0;
  double U               = 2.0;
  double J               = 0.2;
  bool n_n_only          = false;

  // Poles of delta
  std::vector<double> epsilon{2.3, -2.3};

  // Hybridization matrices
  auto V = std::vector{matrix<dcomplex>{num_orbitals, num_orbitals}, matrix<dcomplex>{num_orbitals, num_orbitals}};

  for (int o1 = 0; o1 < num_orbitals; ++o1)
    for (int o2 = 0; o2 < num_orbitals; ++o2) {
      V[0](o1, o2) = (o1 == o2 ? 1.0 : 0.1);
      V[1](o1, o2) = (o1 == o2 ? 1.0 : 0.1);
    }

  assert(epsilon.size() == V.size());

  // Hamiltonian
  many_body_op_t H;
  for (int o = 0; o < num_orbitals; ++o) { H += U * n("up", o) * n("down", o); }
  for (int o1 = 0; o1 < num_orbitals; ++o1)
    for (int o2 = 0; o2 < num_orbitals; ++o2) {
      if (o1 == o2) continue;
      H += (U - 2 * J) * n("up", o1) * n("down", o2);
    }
  for (int o1 = 0; o1 < num_orbitals; ++o1)
    for (int o2 = 0; o2 < num_orbitals; ++o2) {
      if (o2 >= o1) continue;
      H += (U - 3 * J) * n("up", o1) * n("up", o2);
      H += (U - 3 * J) * n("down", o1) * n("down", o2);
    }

  if (!n_n_only) { // spin flips and pair hopping
    for (int o1 = 0; o1 < num_orbitals; ++o1)
      for (int o2 = 0; o2 < num_orbitals; ++o2) {
        if (o1 == o2) continue;
        H += -J * c_dag("up", o1) * c_dag("down", o1) * c("up", o2) * c("down", o2);
        H += -J * c_dag("up", o1) * c_dag("down", o2) * c("up", o2) * c("down", o1);
      }
  }

#ifdef QN
  // quantum numbers
  std::vector<many_body_op_t> qn;
  qn.resize(2);
  for (int o = 0; o < num_orbitals; ++o) {
    qn[0] += n("up", o);
    qn[1] += n("down", o);
  }
#endif

  // gf structure
  indices_t indices{}; 
  for( int o : range(num_orbitals) ) indices.emplace_back(o); 
  gf_struct_t gf_struct{{"up", indices}, {"down", indices}};

  // Construct CTQMC solver
  solver_core solver({beta, gf_struct, 1025, 2500});

  // Set G0
  auto delta_iw = gf<imfreq>{{beta, Fermion}, {num_orbitals, num_orbitals}};

  triqs::clef::placeholder<0> om_;
  auto term = gf<imfreq>{{beta, Fermion}, {num_orbitals, num_orbitals}};
  for (int j = 0; j < epsilon.size(); ++j) {
    term(om_) << 1.0 / (om_ - epsilon[j]);

    matrix<dcomplex> m(num_orbitals, num_orbitals);
    for (std::size_t w_index = 0; w_index < term.mesh().size(); ++w_index) {
      m                                = term.data()(w_index, ellipsis());
      m                                = conj(V[j]) * m * V[j];
      term.data()(w_index, ellipsis()) = m;
    }
    delta_iw = delta_iw + term;
  }

  auto g0_iw = gf<imfreq>{{beta, Fermion}, {num_orbitals, num_orbitals}};
  g0_iw(om_) << om_ + mu - delta_iw(om_);
  solver.G0_iw()[0] = triqs::gfs::inverse(g0_iw);
  solver.G0_iw()[1] = triqs::gfs::inverse(g0_iw);

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
  std::string filename = "kanamori_offdiag";
#ifdef QN
  filename += "_qn";
#endif

  auto &G_tau = *solver.G_tau;

  if (rank == 0) {
    h5::file G_file(filename + ".out.h5", 'w');
    h5_write(G_file, "G_up", G_tau[0]);
    h5_write(G_file, "G_down", G_tau[1]);
  }

  gf<imtime> g;
  if (rank == 0) {
    h5::file G_file(filename + ".ref.h5", 'r');
    h5_read(G_file, "G_up", g);
    EXPECT_GF_NEAR(g, G_tau[0]);
    h5_read(G_file, "G_down", g);
    EXPECT_GF_NEAR(g, G_tau[1]);
  }
}
MAKE_MAIN;
