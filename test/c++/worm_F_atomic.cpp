#include <triqs/test_tools/gfs.hpp>

#include <triqs_cthyb/solver_core.hpp>
#include <triqs_cthyb/impurity_trace.hpp>
#include <triqs_cthyb/measures/F_partition.hpp>

#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <utility>

using namespace triqs_cthyb;
using triqs::operators::c;
using triqs::operators::c_dag;
using triqs::operators::n;
using namespace triqs::gfs;
using namespace triqs::mesh;
using triqs::hilbert_space::fundamental_operator_set;
using triqs::hilbert_space::gf_struct_t;

namespace {

  gf<imtime> exact_atomic_F_tau(triqs_cthyb::atom_diag const &ad, double beta, gf_struct_t const &gf_struct, many_body_op_t const &h_int, int block,
                                int inner_Q, int inner_cdag, int n_tau) {
    auto const &block_name = gf_struct[block].first;

    impurity_trace imp_trace(beta, ad, nullptr);
    auto [atomic_z, tmp] = imp_trace.compute();

    auto Q_op     = h_int * c<h_scalar_t>(block_name, inner_Q) - c<h_scalar_t>(block_name, inner_Q) * h_int;
    auto Q_desc   = imp_trace.attach_aux_operator(Q_op);
    auto cd_desc  = imp_trace.attach_aux_operator(c_dag<h_scalar_t>(block_name, inner_cdag));
    auto tau_seg  = time_segment(beta);
    auto F_tau_ex = gf<imtime>{{beta, Fermion, n_tau}, {1, 1}};

    for (auto tau : F_tau_ex.mesh()) {
      double tau_val = double(tau);
      double eps     = 0.0;
      if (tau_val == 0.0) eps = -1e-14;
      if (tau_val == beta) eps = 1e-14;

      auto tau_cdag = tau_seg.make_time_pt(0.0);
      auto tau_Q    = tau_seg.make_time_pt(tau_val - eps);

      try {
        imp_trace.try_insert(tau_cdag, cd_desc);
        imp_trace.try_insert(tau_Q, Q_desc);
        auto [atomic_weight, atomic_reweighting] = imp_trace.compute();
        F_tau_ex[tau]                            = atomic_weight;
      } catch (rbt_insert_error const &) {
        F_tau_ex[tau] = std::nan("");
      }
      imp_trace.cancel_insert();
    }

    F_tau_ex /= atomic_z;
    return F_tau_ex;
  }

  dcomplex exact_atomic_F_value(triqs_cthyb::atom_diag const &ad, double beta, gf_struct_t const &gf_struct, many_body_op_t const &h_int, int block,
                                int inner_Q, int inner_cdag, double tau_rel, double tau_cdag) {
    auto const &block_name = gf_struct[block].first;

    impurity_trace imp_trace(beta, ad, nullptr);
    auto [atomic_z, tmp] = imp_trace.compute();

    auto Q_op    = h_int * c<h_scalar_t>(block_name, inner_Q) - c<h_scalar_t>(block_name, inner_Q) * h_int;
    auto Q_desc  = imp_trace.attach_aux_operator(Q_op);
    auto cd_desc = imp_trace.attach_aux_operator(c_dag<h_scalar_t>(block_name, inner_cdag));
    auto tau_seg = time_segment(beta);

    double tau_Q = std::fmod(tau_cdag + tau_rel, beta);
    if (tau_Q < 0.0) tau_Q += beta;

    imp_trace.try_insert(tau_seg.make_time_pt(tau_cdag), cd_desc);
    imp_trace.try_insert(tau_seg.make_time_pt(tau_Q), Q_desc);
    auto [atomic_weight, atomic_reweighting] = imp_trace.compute();
    imp_trace.cancel_insert();

    return atomic_weight / atomic_z;
  }

  double hubbard_atom_F_tau(double beta, double U, double mu, double tau) {
    double e_down   = -mu;
    double e_double = -2.0 * mu + U;
    double z        = 1.0 + 2.0 * std::exp(-beta * e_down) + std::exp(-beta * e_double);
    return -U * std::exp(-beta * e_down) * std::exp(tau * (e_down - e_double)) / z;
  }

  struct partition_outputs {
    std::optional<G_tau_t> tau;
    std::optional<G_l_t> legendre;
  };

  partition_outputs run_partition_outputs(bool measure_tau, bool measure_l, long stride) {
    double const beta = 5.0;
    double const U    = 2.0;
    double const mu   = 1.0;
    int const n_iw    = 20;
    int const n_tau   = 41;
    int const n_l     = 8;

    gf_struct_t gf_struct{{"tot", 2}};
    auto h_int = U * n("tot", 0) * n("tot", 1);
    solver_core solver({.beta = beta, .gf_struct = gf_struct, .n_iw = n_iw, .n_tau = n_tau, .n_l = n_l});

    nda::clef::placeholder<0> om_;
    auto delta_iw = gf<imfreq>{{beta, Fermion, n_iw}, {2, 2}};
    nda::matrix<dcomplex> bath_coupling(2, 2);
    bath_coupling()  = 0.0;
    bath_coupling(0, 0) = 1.0;
    bath_coupling(1, 1) = 0.7;
    delta_iw(om_) << bath_coupling * (1.0 / (om_ - 2.0) + 1.0 / (om_ + 2.0));

    auto g0_iw = gf<imfreq>{{beta, Fermion, n_iw}, {2, 2}};
    g0_iw(om_) << om_ + mu - delta_iw(om_);
    solver.G0_iw()[0] = triqs::gfs::inverse(g0_iw);

    auto p                       = solve_parameters_t{.h_int = h_int, .n_cycles = 1200};
    p.length_cycle               = 5;
    p.n_warmup_cycles            = 100;
    p.random_seed                = 24680;
    p.random_name                = "";
    p.verbosity                  = 0;
    p.move_double                = false;
    p.partition_method           = "none";
    p.measure_G_tau              = false;
    p.measure_F_tau_partition    = measure_tau;
    p.measure_F_l_partition      = measure_l;
    p.measure_F_partition_stride = stride;
    solver.solve(p);

    return {solver.F_tau_partition, solver.F_l_partition};
  }

} // namespace

TEST(PartitionF, StrideScheduleAndHdfRoundTrip) {
  triqs_cthyb::detail::partition_measurement_schedule schedule(3);
  for (int event = 0; event < 10; ++event) EXPECT_EQ(schedule.select_next(), event % 3 == 0);
  EXPECT_ANY_THROW(triqs_cthyb::detail::partition_measurement_schedule(0));
  EXPECT_ANY_THROW(triqs_cthyb::detail::partition_measurement_schedule(-1));

  solve_parameters_t parameters;
  EXPECT_EQ(parameters.measure_F_partition_stride, 1);
  parameters.measure_F_partition_stride = 7;
  auto const filename = "partition_stride." + std::to_string(mpi::communicator().rank()) + ".h5";
  {
    h5::file file(filename, 'w');
    h5_write(file, "parameters", parameters);
  }
  solve_parameters_t restored;
  {
    h5::file file(filename, 'r');
    h5_read(file, "parameters", restored);
  }
  std::remove(filename.c_str());
  EXPECT_EQ(restored.measure_F_partition_stride, 7);
}

TEST(PartitionF, FusedOptionalOutputsMatchSingleOutputMeasures) {
  auto tau_only = run_partition_outputs(true, false, 3);
  auto l_only   = run_partition_outputs(false, true, 3);
  auto both     = run_partition_outputs(true, true, 3);

  ASSERT_TRUE(tau_only.tau.has_value());
  EXPECT_FALSE(tau_only.legendre.has_value());
  EXPECT_FALSE(l_only.tau.has_value());
  ASSERT_TRUE(l_only.legendre.has_value());
  ASSERT_TRUE(both.tau.has_value());
  ASSERT_TRUE(both.legendre.has_value());

  EXPECT_BLOCK_GF_NEAR(*both.tau, *tau_only.tau, 1e-14);
  EXPECT_BLOCK_GF_NEAR(*both.legendre, *l_only.legendre, 1e-14);
}

TEST(WormF, HubbardAtomExactTrace) {
  double beta = 5.0;
  double U    = 2.0;
  double mu   = 1.0;
  int n_tau   = 201;

  gf_struct_t gf_struct{{"up", 1}, {"down", 1}};
  auto fops  = fundamental_operator_set(gf_struct);
  auto h_int = U * n("up", 0) * n("down", 0);
  auto h_loc = h_int - mu * (n("up", 0) + n("down", 0));
  auto ad    = triqs_cthyb::atom_diag(h_loc, fops);

  auto F_up = exact_atomic_F_tau(ad, beta, gf_struct, h_int, 0, 0, 0, n_tau);

  for (auto tau : F_up.mesh()) {
    double tau_val = double(tau);
    EXPECT_NEAR(real(F_up[tau](0, 0)), hubbard_atom_F_tau(beta, U, mu, tau_val), 1e-12);
  }
}

TEST(WormF, HubbardAtomAbsoluteTimeTraceSign) {
  double beta = 5.0;
  double U    = 2.0;
  double mu   = 1.0;

  gf_struct_t gf_struct{{"up", 1}, {"down", 1}};
  auto fops  = fundamental_operator_set(gf_struct);
  auto h_int = U * n("up", 0) * n("down", 0);
  auto h_loc = h_int - mu * (n("up", 0) + n("down", 0));
  auto ad    = triqs_cthyb::atom_diag(h_loc, fops);

  for (auto [tau_rel, tau_cdag] : {std::pair{0.25, 0.10}, std::pair{0.25, 4.90}, std::pair{2.50, 0.75}, std::pair{4.25, 2.00}}) {
    auto F = exact_atomic_F_value(ad, beta, gf_struct, h_int, 0, 0, 0, tau_rel, tau_cdag);
    EXPECT_NEAR(real(F), hubbard_atom_F_tau(beta, U, mu, tau_rel), 1e-12);
  }
}

TEST(WormF, HubbardAtomStochasticEstimatorSmoke) {
  double beta = 5.0;
  double U    = 2.0;
  double mu   = 1.0;
  int n_iw    = 80;
  int n_tau   = 201;

  gf_struct_t gf_struct{{"up", 1}, {"down", 1}};
  auto h_int = U * n("up", 0) * n("down", 0);

  solver_core solver({.beta = beta, .gf_struct = gf_struct, .n_iw = n_iw, .n_tau = n_tau});

  nda::clef::placeholder<0> om_;
  auto g0_iw = gf<imfreq>{{beta, Fermion, n_iw}, {1, 1}};
  g0_iw(om_) << om_ + mu;
  solver.G0_iw()[0] = triqs::gfs::inverse(g0_iw);
  solver.G0_iw()[1] = triqs::gfs::inverse(g0_iw);

  auto p              = solve_parameters_t{.h_int = h_int, .n_cycles = 300000};
  p.length_cycle      = 20;
  p.n_warmup_cycles   = 10000;
  p.random_seed       = 111;
  p.random_name       = "";
  p.verbosity         = 0;
  p.move_double       = false;
  p.use_norm_as_weight = true;
  p.measure_G_tau     = false;
  p.measure_density_matrix = true;
  p.measure_F_tau_worm = true;
  p.worm_eta          = 0.1;
  p.worm_prob         = 10.0;

  solver.solve(p);

  EXPECT_NEAR(real(solver.average_sign()), 1.0, 1e-12);
  EXPECT_EQ(solver.average_sign(), solver.average_sign_partition());
  EXPECT_NEAR(real(solver.average_sign_worm()), -1.0, 1e-12);

  EXPECT_FALSE(solver.density_matrix().empty());
  ASSERT_TRUE(solver.F_tau.has_value());
  auto const &F_up = (*solver.F_tau)[0];

  auto sample = [&F_up](double tau) { return real(F_up[closest_mesh_pt(tau)](0, 0)); };

  auto f0   = sample(0.0);
  auto f025 = sample(0.25);
  auto f05  = sample(0.50);
  auto f125 = sample(1.25);
  auto f25  = sample(2.50);
  auto fb   = sample(beta);

  EXPECT_NEAR(f0 + fb, -1.0, 0.40);
  EXPECT_NEAR(f025, hubbard_atom_F_tau(beta, U, mu, 0.25), 0.25);
  EXPECT_NEAR(f05, hubbard_atom_F_tau(beta, U, mu, 0.50), 0.20);
  EXPECT_NEAR(f125, hubbard_atom_F_tau(beta, U, mu, 1.25), 0.18);
  EXPECT_LT(f025, f05);
  EXPECT_LT(f05, f125);
  EXPECT_LT(f125, f25);
  EXPECT_LT(f25, 0.0);
}

TEST(PartitionF, NormReweightingZeroTraceConfigurationsStayFinite) {
  double beta = 5.0;
  double U    = 2.0;
  double mu   = 1.0;
  int n_iw    = 40;
  int n_tau   = 101;
  int n_l     = 20;

  gf_struct_t gf_struct{{"tot", 2}};
  auto h_int = U * n("tot", 0) * n("tot", 1);

  solver_core solver({.beta = beta, .gf_struct = gf_struct, .n_iw = n_iw, .n_tau = n_tau, .n_l = n_l});

  nda::clef::placeholder<0> om_;
  auto delta_iw = gf<imfreq>{{beta, Fermion, n_iw}, {2, 2}};
  nda::matrix<dcomplex> bath_coupling(2, 2);
  bath_coupling(0, 0) = 1.0;
  bath_coupling(0, 1) = 1.0;
  bath_coupling(1, 0) = 1.0;
  bath_coupling(1, 1) = 1.0;
  delta_iw(om_) << bath_coupling * (1.0 / (om_ - 2.0) + 1.0 / (om_ + 2.0));

  auto g0_iw = gf<imfreq>{{beta, Fermion, n_iw}, {2, 2}};
  g0_iw(om_) << om_ + mu - delta_iw(om_);
  solver.G0_iw()[0] = triqs::gfs::inverse(g0_iw);

  auto p                    = solve_parameters_t{.h_int = h_int, .n_cycles = 50000};
  p.length_cycle            = 10;
  p.n_warmup_cycles         = 5000;
  p.random_seed             = 8675309;
  p.random_name             = "";
  p.verbosity               = 0;
  p.move_double             = false;
  p.partition_method        = "none";
  p.use_norm_as_weight      = true;
  p.measure_G_tau           = false;
  p.measure_F_tau_partition = true;
  p.measure_F_l_partition   = true;
  p.measure_F_partition_stride = 3;

  solver.solve(p);

  auto count_nonfinite = [](auto const &block_gf) {
    long count = 0;
    for (auto const &block : block_gf)
      for (auto const &value : block.data()) count += !isfinite(value);
    return count;
  };

  ASSERT_TRUE(solver.F_tau_partition.has_value());
  EXPECT_EQ(count_nonfinite(*solver.F_tau_partition), 0);

  ASSERT_TRUE(solver.F_l_partition.has_value());
  EXPECT_EQ(count_nonfinite(*solver.F_l_partition), 0);
}

MAKE_MAIN;
