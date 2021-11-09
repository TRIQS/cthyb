// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

#include <triqs/test_tools/gfs.hpp>

#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp> // gf_struct_t
using gf_struct_t = triqs::hilbert_space::gf_struct_t;

using namespace nda;
using namespace triqs::hilbert_space;
using namespace triqs::atom_diag;
using namespace triqs::operators;

// -----------------------------------------------------------------------------

#include <triqs_cthyb/types.hpp>
#include <triqs_cthyb/impurity_trace.hpp>
#include <triqs_cthyb/configuration.hpp> // for op_desc

using linindex_t = std::map<std::pair<int, int>, int>;

// -----------------------------------------------------------------------------
linindex_t make_linear_index(const gf_struct_t &gf_struct, const fundamental_operator_set &fops) {
  linindex_t linindex;
  int block_index = 0;
  for (auto const &[bl, bl_size] : gf_struct) {
    for (auto a : range(bl_size)) { linindex[std::make_pair(block_index, a)] = fops[{bl, a}]; }
    block_index++;
  }
  return linindex;
}

// -----------------------------------------------------------------------------
TEST(impurity_trace, try_insert_cancel_bug) {

  gf_struct_t gf_struct{{"up", 1}, {"dn", 1}};
  fundamental_operator_set fops(gf_struct);
  auto linindex  = make_linear_index(gf_struct, fops);

  // -----------------------------------------------------------------------------
  // atom_diag

  double U  = 1.0;
  double mu = 0.5 * U;

  many_body_operator_real H;
  H +=  - mu * ( n("up", 0 ) + n("dn", 0) ) + U * n("up", 0 ) * n("dn", 0);

  auto ad   = triqs::atom_diag::atom_diag<triqs_cthyb::is_h_scalar_complex>(H, fops);
  std::cout << "Found " << ad.n_subspaces() << " subspaces." << std::endl;

  // -----------------------------------------------------------------------------
  // impurity_trace

  double beta = 1.0;
  triqs_cthyb::impurity_trace imp_trace(beta, ad, nullptr);

  // -----------------------------------------------------------------------------
  // test insertion

  int oidx = 0;
  int block_index = 0;
  
  auto op1        = triqs_cthyb::op_desc{block_index, oidx, true, linindex[std::make_pair(block_index, oidx)]};
  auto op2        = triqs_cthyb::op_desc{block_index, oidx, false, linindex[std::make_pair(block_index, oidx)]};

  triqs_cthyb::time_segment tau_seg(beta);

  // Generate insert error by inserting at the same time  
  auto tau1 = tau_seg.make_time_pt(0.);
  auto tau2 = tau_seg.make_time_pt(0.);

  try {
    
    std::cout << "Before insert: \n" << imp_trace << "\n";

    imp_trace.try_insert(tau1, op1);

    std::cout << "After 1st insert: \n" << imp_trace << "\n";

    imp_trace.try_insert(tau2, op2);

    std::cout << "After 2nd insert: \n" << imp_trace << "\n";

  } catch (rbt_insert_error const &) {

    std::cerr << "Insert error : recovering ... " << std::endl;

    std::cout << "Before cancel: \n" << imp_trace << "\n";
    assert( imp_trace.tree_size == 1 );

    imp_trace.cancel_insert();

    std::cout << "After cancel: \n" << imp_trace << "\n";
    assert( imp_trace.tree_size == 0 );
  }
}

MAKE_MAIN;
