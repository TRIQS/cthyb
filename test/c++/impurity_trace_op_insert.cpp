// -----------------------------------------------------------------------------

#include <triqs/test_tools/gfs.hpp>

#include <triqs/atom_diag/atom_diag.hpp>
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

using h_scalar_t = double;
using linindex_t = std::map<std::pair<int, int>, int>;

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
TEST(atom_diag, op_matrix) {

  gf_struct_t gf_struct{{"up", {0}}, {"dn", {0}}};
  fundamental_operator_set fops(gf_struct);
  auto linindex = make_linear_index(gf_struct, fops);

  // -----------------------------------------------------------------------------
  // atom_diag

  double U  = 1.0;
  double mu = 0.1 * U;

  many_body_operator_real H;
  H += -mu * (n("up", 0) + n("dn", 0)) + U * n("up", 0) * n("dn", 0);

  auto ad = triqs::atom_diag::atom_diag<triqs_cthyb::is_h_scalar_complex>(H, fops);
  std::cout << "Found " << ad.n_subspaces() << " subspaces." << std::endl;

  // -----------------------------------------------------------------------------

  // -----------------------------------------------------------------------------
  // impurity_trace

  double beta = 2.0;
  triqs_cthyb::impurity_trace imp_trace(beta, ad, nullptr);

  triqs_cthyb::h_scalar_t atomic_z, tmp;
  std::tie(atomic_z, tmp) = imp_trace.compute();

  std::cout << "Z = " << atomic_z << "\n";

  // -----------------------------------------------------------------------------

  triqs_cthyb::time_segment tau_seg(beta);
  triqs_cthyb::h_scalar_t new_atomic_weight, new_atomic_reweighting;

  // -----------------------------------------------------------------------------

  {
    many_body_operator_real op = n("dn", 0) * n("up", 0);
    auto op_d                  = imp_trace.attach_aux_operator(op);
    auto tau1                  = tau_seg.make_time_pt(0.);

    try {
      imp_trace.try_insert(tau1, op_d);
      std::cout << imp_trace << "\n";
      std::tie(new_atomic_weight, new_atomic_reweighting) = imp_trace.compute();
    } catch (rbt_insert_error const &) {
      std::cerr << "Insert error : recovering ... " << std::endl;
      new_atomic_weight      = std::nan("");
      new_atomic_reweighting = std::nan("");
    }

    imp_trace.cancel_insert();
    std::cout << new_atomic_weight << ", " << new_atomic_reweighting << "\n";

    triqs_cthyb::h_scalar_t exp_val = new_atomic_weight / atomic_z;
    std::cout << "exp_val = " << exp_val << "\n";
  }

  // -----------------------------------------------------------------------------
  // gf eval, using the imp_trace

  int ntau = 10;
  auto g = gf<imtime>{{beta, Fermion, ntau}, {1, 1}};

  many_body_operator_real op1 = c_dag("up", 0);
  many_body_operator_real op2 = n("dn", 0) * c("up", 0);

  auto op1_d = imp_trace.attach_aux_operator(op1);
  auto op2_d = imp_trace.attach_aux_operator(op2);

  for (auto tau : g.mesh()) {

    double eps = 0;
    if (tau == 0. ) eps = -1e-14; // This should not be needed FIXME
    if (tau == beta) eps = 1e-14; // This should not be needed FIXME

    auto tau1 = tau_seg.make_time_pt(0.);
    auto tau2 = tau_seg.make_time_pt(tau - eps);
    
    try {
      imp_trace.try_insert(tau1, op1_d);
      imp_trace.try_insert(tau2, op2_d);
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
  
  {
    h5::file fd("impurity_trace_op_insert.h5", 'w');
    h5_write(fd, "g", g);
  }  
  
}

MAKE_MAIN;
