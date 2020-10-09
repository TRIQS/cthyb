// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

#include <cmath>

#include <triqs/test_tools/gfs.hpp>

#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/atom_diag/functions.hpp>
#include <triqs/atom_diag/gf.hpp>
#include <h5/h5.hpp>

#include <triqs/hilbert_space/fundamental_operator_set.hpp> // gf_struct_t
using gf_struct_t = triqs::hilbert_space::gf_struct_t;

using namespace triqs::arrays;
using namespace triqs::hilbert_space;
using namespace triqs::atom_diag;
using namespace triqs::operators;

// -----------------------------------------------------------------------------

#include <triqs_cthyb/types.hpp>
#include <triqs_cthyb/impurity_trace.hpp>
#include <triqs_cthyb/configuration.hpp> // for op_desc

using linindex_t = std::map<std::pair<int, int>, int>;

// -----------------------------------------------------------------------------
gf_struct_t make_gf_struct(int n_orb) {
  indices_t indices{};
  for (int o : range(n_orb)) indices.emplace_back(o);
  gf_struct_t gf_struct{{"up", indices}, {"dn", indices}};
  return gf_struct;
}

// -----------------------------------------------------------------------------
linindex_t make_linear_index(const gf_struct_t &gf_struct, const fundamental_operator_set &fops) {
  linindex_t linindex;
  int block_index = 0;
  for (auto const &bl : gf_struct) {
    int inner_index = 0;
    for (auto const &a : bl.second) {
      linindex[std::make_pair(block_index, inner_index)] = fops[{bl.first, a}];
      inner_index++;
    }
    block_index++;
  }
  return linindex;
}

// -----------------------------------------------------------------------------
template <typename OP> OP make_hamiltonian(int n_orb, double mu, double U, double J) {

  auto orbs = range(n_orb);

  OP h;

  for (int o : orbs) h += -mu * (n("up", o) + n("dn", o));

  // Density-density interactions
  for (int o : orbs) h += U * n("up", o) * n("dn", o);

  for (int o1 : orbs) {
    for (int o2 : orbs) {
      if (o1 == o2) continue;
      h += (U - 2 * J) * n("up", o1) * n("dn", o2);
    }
  }

  for (int o1 : orbs)
    for (int o2 : orbs) {
      if (o2 >= o1) continue;
      h += (U - 3 * J) * n("up", o1) * n("up", o2);
      h += (U - 3 * J) * n("dn", o1) * n("dn", o2);
    }

  // spin-flip and pair-hopping
  for (int o1 : orbs) {
    for (int o2 : orbs) {
      if (o1 == o2) continue;
      h += -J * c_dag("up", o1) * c_dag("dn", o1) * c("up", o2) * c("dn", o2);
      h += -J * c_dag("up", o1) * c_dag("dn", o2) * c("up", o2) * c("dn", o1);
    }
  }

  return h;
}

// -----------------------------------------------------------------------------
TEST(impurity_trace, atomic_gf) {

  int n_orb = 3;

  auto gf_struct = make_gf_struct(n_orb);
  auto fops      = fundamental_operator_set(gf_struct);
  auto linindex  = make_linear_index(gf_struct, fops);

  // -----------------------------------------------------------------------------
  // atom_diag

  double U  = 1.0;
  double J  = 0.2;
  double mu = 0.5 * U;
  auto H    = make_hamiltonian<many_body_operator_real>(n_orb, mu, U, J);
  auto ad   = triqs::atom_diag::atom_diag<triqs_cthyb::is_h_scalar_complex>(H, fops);
  std::cout << "Found " << ad.n_subspaces() << " subspaces." << std::endl;

  // -----------------------------------------------------------------------------
  // impurity_trace

  double beta = 1.0;
  triqs_cthyb::impurity_trace imp_trace(beta, ad, nullptr);

  triqs_cthyb::h_scalar_t atomic_z, tmp;
  std::tie(atomic_z, tmp) = imp_trace.compute();

  // -----------------------------------------------------------------------------
  // test insertion

  int block_index = 0;
  int oidx        = 0;
  auto op1        = triqs_cthyb::op_desc{block_index, oidx, true, linindex[std::make_pair(block_index, oidx)]};
  auto op2        = triqs_cthyb::op_desc{block_index, oidx, false, linindex[std::make_pair(block_index, oidx)]};

  triqs_cthyb::time_segment tau_seg(beta);

  triqs_cthyb::h_scalar_t new_atomic_weight, new_atomic_reweighting;

  // -----------------------------------------------------------------------------
  // gf eval, using the imp_trace

  int ntau = 10;

  auto g = gf<imtime>{{beta, Fermion, ntau}, {1, 1}};

  for (auto tau : g.mesh()) {

    double eps = 0;
    if (tau == 0. ) eps = -1e-14; // This should not be needed FIXME
    if (tau == beta) eps = 1e-14; // This should not be needed FIXME

    auto tau1 = tau_seg.make_time_pt(0.);
    auto tau2 = tau_seg.make_time_pt(tau - eps);
    
    try {
      imp_trace.try_insert(tau1, op1);
      imp_trace.try_insert(tau2, op2);
      std::tie(new_atomic_weight, new_atomic_reweighting) = imp_trace.compute();
    } catch (rbt_insert_error const &) {
      std::cerr << "Insert error : recovering ... " << std::endl;
      new_atomic_weight = std::nan("");
      new_atomic_reweighting = std::nan("");
    }
    
    imp_trace.cancel_insert();

    g[tau] = new_atomic_weight;
  }

  g /= -atomic_z;

  // -----------------------------------------------------------------------------
  // gf eval, using atom_diag
  
  auto g_ref = atomic_g_tau(ad, beta, gf_struct, ntau, {});

  // compare single particle Green's function from atom_diag with result from trace
  EXPECT_GF_NEAR(g, slice_target(g_ref[0], range(0,1), range(0,1)));
  
  {
    h5::file fd("impurity_trace_atomic_gf.h5", 'w');
    h5_write(fd, "g", g);
    h5_write(fd, "g_ref", g_ref);
  }
}

MAKE_MAIN;
