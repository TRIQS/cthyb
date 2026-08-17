// -----------------------------------------------------------------------------

#include <cmath>

#include <triqs/test_tools/gfs.hpp>
#include <triqs/atom_diag/atom_diag.hpp>
#include <triqs/gfs/block/gf_struct.hpp>

using namespace nda;
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
  for (auto const &[bl, bl_size] : gf_struct) {
    for (auto a : range(bl_size)) { linindex[std::make_pair(block_index, a)] = fops[{bl, a}]; }
    block_index++;
  }
  return linindex;
}

// -----------------------------------------------------------------------------
TEST(atom_diag, op_matrix) {

  gf_struct_t gf_struct{{"up", 1}, {"dn", 1}};
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

TEST(impurity_trace, delete_by_time_key) {

  gf_struct_t gf_struct{{"up", 1}, {"dn", 1}};
  fundamental_operator_set fops(gf_struct);

  double U  = 1.0;
  double mu = 0.1 * U;

  many_body_operator_real H;
  H += -mu * (n("up", 0) + n("dn", 0)) + U * n("up", 0) * n("dn", 0);

  auto ad = triqs::atom_diag::atom_diag<triqs_cthyb::is_h_scalar_complex>(H, fops);

  double beta = 2.0;
  triqs_cthyb::time_segment tau_seg(beta);
  auto tau_up = tau_seg.make_time_pt(0.4);
  auto tau_dn = tau_seg.make_time_pt(1.3);

  triqs_cthyb::impurity_trace imp_trace(beta, ad, nullptr);
  auto n_up = imp_trace.attach_aux_operator(n("up", 0));
  auto n_dn = imp_trace.attach_aux_operator(n("dn", 0));
  EXPECT_EQ(n_up.block_index, -1);
  EXPECT_EQ(n_dn.block_index, -1);

  imp_trace.try_insert(tau_up, n_up);
  imp_trace.try_insert(tau_dn, n_dn);
  imp_trace.confirm_insert();
  auto [both_weight, both_reweighting] = imp_trace.compute();

  triqs_cthyb::impurity_trace ref_both(beta, ad, nullptr);
  auto ref_n_up = ref_both.attach_aux_operator(n("up", 0));
  auto ref_n_dn = ref_both.attach_aux_operator(n("dn", 0));
  ref_both.try_insert(tau_up, ref_n_up);
  ref_both.try_insert(tau_dn, ref_n_dn);
  ref_both.confirm_insert();
  auto [ref_both_weight, ref_both_reweighting] = ref_both.compute();

  EXPECT_LT(std::abs(both_weight - ref_both_weight), 1e-12);
  EXPECT_LT(std::abs(both_reweighting - ref_both_reweighting), 1e-12);

  triqs_cthyb::impurity_trace ref_deleted(beta, ad, nullptr);
  auto ref_deleted_n_dn = ref_deleted.attach_aux_operator(n("dn", 0));
  ref_deleted.try_insert(tau_dn, ref_deleted_n_dn);
  ref_deleted.confirm_insert();
  auto [ref_deleted_weight, ref_deleted_reweighting] = ref_deleted.compute();

  imp_trace.try_delete(tau_up);
  auto [trial_weight, trial_reweighting] = imp_trace.compute();
  EXPECT_LT(std::abs(trial_weight - ref_deleted_weight), 1e-12);
  EXPECT_LT(std::abs(trial_reweighting - ref_deleted_reweighting), 1e-12);

  imp_trace.cancel_delete();
  auto [cancel_weight, cancel_reweighting] = imp_trace.compute();
  EXPECT_LT(std::abs(cancel_weight - ref_both_weight), 1e-12);
  EXPECT_LT(std::abs(cancel_reweighting - ref_both_reweighting), 1e-12);

  imp_trace.try_delete(tau_up);
  imp_trace.confirm_delete();
  auto [deleted_weight, deleted_reweighting] = imp_trace.compute();
  EXPECT_LT(std::abs(deleted_weight - ref_deleted_weight), 1e-12);
  EXPECT_LT(std::abs(deleted_reweighting - ref_deleted_reweighting), 1e-12);
}

TEST(impurity_trace, single_key_replace_cancel_confirm) {

  gf_struct_t gf_struct{{"up", 1}, {"dn", 1}};
  fundamental_operator_set fops(gf_struct);
  auto linindex = make_linear_index(gf_struct, fops);

  many_body_operator_real H;
  H += -0.2 * (n("up", 0) + n("dn", 0)) + 1.3 * n("up", 0) * n("dn", 0);
  H += 0.17 * (c_dag("up", 0) * c("dn", 0) + c_dag("dn", 0) * c("up", 0));
  auto ad = triqs::atom_diag::atom_diag<triqs_cthyb::is_h_scalar_complex>(H, fops);

  double beta = 2.0;
  triqs_cthyb::time_segment tau_seg(beta);
  auto tau_0 = tau_seg.make_time_pt(0.15);
  auto tau_1 = tau_seg.make_time_pt(0.45);
  auto tau_2 = tau_seg.make_time_pt(0.9);
  auto tau_3 = tau_seg.make_time_pt(1.55);

  triqs_cthyb::op_desc cdag_up{0, 0, true, linindex.at(std::make_pair(0, 0))};
  triqs_cthyb::op_desc c_up{0, 0, false, linindex.at(std::make_pair(0, 0))};
  triqs_cthyb::op_desc cdag_dn{1, 0, true, linindex.at(std::make_pair(1, 0))};
  triqs_cthyb::op_desc c_dn{1, 0, false, linindex.at(std::make_pair(1, 0))};

  triqs_cthyb::impurity_trace imp_trace(beta, ad, nullptr, true);
  imp_trace.try_insert(tau_0, cdag_up);
  imp_trace.try_insert(tau_1, c_up);
  imp_trace.try_insert(tau_2, cdag_dn);
  imp_trace.try_insert(tau_3, c_dn);
  imp_trace.confirm_insert();

  auto scaled_c = imp_trace.attach_aux_operator(2.5 * c("up", 0));
  auto physical_trace = [&]() {
    auto [weight, reweighting] = imp_trace.compute();
    return weight * reweighting;
  };

  auto baseline = physical_trace();
  triqs_cthyb::configuration::oplist_t generic_update;
  generic_update.emplace(tau_1, scaled_c);
  imp_trace.try_replace(generic_update);
  auto generic_trial = physical_trace();
  imp_trace.cancel_replace();
  EXPECT_LE(std::abs(physical_trace() - baseline), 1.e-11 * std::max(1.0, std::abs(baseline)));

  imp_trace.try_replace(tau_1, scaled_c);
  auto single_trial = physical_trace();
  EXPECT_LE(std::abs(single_trial - generic_trial), 1.e-11 * std::max(1.0, std::abs(generic_trial)));
  EXPECT_ANY_THROW(imp_trace.try_replace(tau_2, cdag_dn));
  imp_trace.cancel_replace();
  EXPECT_LE(std::abs(physical_trace() - baseline), 1.e-11 * std::max(1.0, std::abs(baseline)));

  auto missing = tau_seg.make_time_pt(1.1);
  EXPECT_ANY_THROW(imp_trace.try_replace(missing, scaled_c));
  EXPECT_LE(std::abs(physical_trace() - baseline), 1.e-11 * std::max(1.0, std::abs(baseline)));

  imp_trace.try_replace(tau_1, scaled_c);
  imp_trace.confirm_replace();
  EXPECT_LE(std::abs(physical_trace() - generic_trial), 1.e-11 * std::max(1.0, std::abs(generic_trial)));
  imp_trace.cancel_replace(); // no active replacement: must be a safe no-op
  EXPECT_LE(std::abs(physical_trace() - generic_trial), 1.e-11 * std::max(1.0, std::abs(generic_trial)));

  triqs_cthyb::configuration::oplist_t multi_update;
  multi_update.emplace(tau_0, cdag_dn);
  multi_update.emplace(tau_3, c_up);
  imp_trace.try_replace(multi_update);
  auto multi_trial = physical_trace();
  imp_trace.cancel_replace();
  EXPECT_LE(std::abs(physical_trace() - generic_trial), 1.e-11 * std::max(1.0, std::abs(generic_trial)));
  imp_trace.try_replace(multi_update);
  imp_trace.confirm_replace();
  EXPECT_LE(std::abs(physical_trace() - multi_trial), 1.e-11 * std::max(1.0, std::abs(multi_trial)));
}

TEST(impurity_trace, batch_single_replacements_match_generic_reference) {

  gf_struct_t gf_struct{{"up", 1}, {"dn", 1}};
  fundamental_operator_set fops(gf_struct);
  auto linindex = make_linear_index(gf_struct, fops);

  many_body_operator_real H;
  H += -0.2 * (n("up", 0) + n("dn", 0)) + 1.3 * n("up", 0) * n("dn", 0);
  H += 0.17 * (c_dag("up", 0) * c("dn", 0) + c_dag("dn", 0) * c("up", 0));
  auto ad = triqs::atom_diag::atom_diag<triqs_cthyb::is_h_scalar_complex>(H, fops);

  double beta = 2.0;
  triqs_cthyb::time_segment tau_seg(beta);
  std::vector<triqs_cthyb::time_pt> times;
  for (double tau : {0.12, 0.36, 0.62, 0.88, 1.18, 1.61}) times.push_back(tau_seg.make_time_pt(tau));

  triqs_cthyb::op_desc cdag_up{0, 0, true, linindex.at(std::make_pair(0, 0))};
  triqs_cthyb::op_desc c_up{0, 0, false, linindex.at(std::make_pair(0, 0))};
  triqs_cthyb::op_desc cdag_dn{1, 0, true, linindex.at(std::make_pair(1, 0))};
  triqs_cthyb::op_desc c_dn{1, 0, false, linindex.at(std::make_pair(1, 0))};

  for (bool use_norm_as_weight : {false, true}) {
      triqs_cthyb::impurity_trace imp_trace(beta, ad, nullptr, use_norm_as_weight, true);
      for (auto const &[tau, op] : std::vector<std::pair<triqs_cthyb::time_pt, triqs_cthyb::op_desc>>{
              {times[0], cdag_up}, {times[1], c_up}, {times[2], cdag_dn}, {times[3], c_dn}, {times[4], cdag_up}, {times[5], c_up}}) {
        imp_trace.try_insert(tau, op);
        imp_trace.confirm_insert();
      }

      auto proportional = imp_trace.attach_aux_operator(2.5 * c("up", 0));
      auto general      = imp_trace.attach_aux_operator(n("dn", 0) * c("up", 0));
      auto different    = imp_trace.attach_aux_operator(c("up", 0));

      using request_t = triqs_cthyb::impurity_trace::single_replacement_request;
      std::vector<request_t> requests{{times[5], general}, {times[1], proportional}, {times[3], different}};

      auto [baseline_weight, baseline_reweighting] = imp_trace.compute();
      auto baseline_trace                          = baseline_weight * baseline_reweighting;
      auto batch                                   = imp_trace.compute_single_replacement_traces(requests);
      auto repeated_batch                          = imp_trace.compute_single_replacement_traces(requests);
      ASSERT_EQ(batch.size(), requests.size());
      ASSERT_EQ(repeated_batch.size(), requests.size());

      for (std::size_t i = 0; i < requests.size(); ++i) {
        triqs_cthyb::configuration::oplist_t generic_update;
        generic_update.emplace(requests[i].key, requests[i].replacement);
        imp_trace.try_replace(generic_update);
        auto [reference_weight, reference_reweighting] = imp_trace.compute();
        auto reference_trace                           = reference_weight * reference_reweighting;
        imp_trace.cancel_replace();

        auto tolerance = 1.e-11 * std::max(1.0, std::abs(reference_trace));
        EXPECT_LE(std::abs(batch[i] - reference_trace), tolerance);
        EXPECT_LE(std::abs(repeated_batch[i] - reference_trace), tolerance);
      }

      auto counts = imp_trace.get_last_single_replacement_path_counts();
      EXPECT_GT(counts.structural_zero, 0);
      EXPECT_GT(counts.proportional, 0);
      EXPECT_GT(counts.general, 0);

      auto [after_weight, after_reweighting] = imp_trace.compute();
      EXPECT_LE(std::abs(after_weight * after_reweighting - baseline_trace), 1.e-11 * std::max(1.0, std::abs(baseline_trace)));

      std::vector<request_t> duplicate{{times[1], proportional}, {times[1], general}};
      EXPECT_ANY_THROW((void)imp_trace.compute_single_replacement_traces(duplicate));
      std::vector<request_t> missing{{tau_seg.make_time_pt(1.02), proportional}};
      EXPECT_ANY_THROW((void)imp_trace.compute_single_replacement_traces(missing));
      std::vector<request_t> invalid{{times[1], triqs_cthyb::op_desc{-1, 0, true, -1000}}};
      EXPECT_ANY_THROW((void)imp_trace.compute_single_replacement_traces(invalid));

      imp_trace.try_replace(times[1], proportional);
      EXPECT_ANY_THROW((void)imp_trace.compute_single_replacement_traces(requests));
      imp_trace.cancel_replace();
      imp_trace.try_delete(times[3]);
      EXPECT_ANY_THROW((void)imp_trace.compute_single_replacement_traces(requests));
      imp_trace.cancel_delete();
      auto pending_tau = tau_seg.make_time_pt(1.91);
      imp_trace.try_insert(pending_tau, cdag_up);
      EXPECT_ANY_THROW((void)imp_trace.compute_single_replacement_traces(requests));
      imp_trace.cancel_insert();
      auto [restored_weight, restored_reweighting] = imp_trace.compute();
      EXPECT_LE(std::abs(restored_weight * restored_reweighting - baseline_trace), 1.e-11 * std::max(1.0, std::abs(baseline_trace)));
  }
}

MAKE_MAIN;
