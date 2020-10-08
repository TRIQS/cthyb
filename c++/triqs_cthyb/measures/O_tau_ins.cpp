/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2018, The Simons Foundation
 * Author: H. U.R. Strand
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <triqs/mc_tools.hpp>

#include "./O_tau_ins.hpp"

namespace triqs_cthyb {

  using namespace triqs::gfs;
  using namespace triqs::mesh;

  measure_O_tau_ins::measure_O_tau_ins(std::optional<gf<imtime, scalar_valued>> &O_tau_opt, qmc_data const &data, int n_tau,
                                       many_body_op_t const &op1, many_body_op_t const &op2, int min_ins, mc_tools::random_generator &rng)
    : data(data), average_sign(0), op1(op1), op2(op2), min_ins(min_ins), rng(rng) {
    O_tau_opt = gf<imtime, scalar_valued>{{data.config.beta(), Boson, n_tau}};
    O_tau.rebind(*O_tau_opt);
    O_tau() = 0.0;

    op1_d = data.imp_trace.attach_aux_operator(op1);
    op2_d = data.imp_trace.attach_aux_operator(op2);
  }

  void measure_O_tau_ins::accumulate(mc_weight_t s) {
    s *= data.atomic_reweighting;
    average_sign += s;

    int pto = 0;
    for (const auto &det : data.dets) pto += det.size();
    int nsamples = pto * pto;
    if( nsamples < min_ins ) nsamples = min_ins;

    mc_weight_t atomic_weight, atomic_reweighting;
    auto [bare_atomic_weight, bare_atomic_reweighting] = data.imp_trace.compute();
    const auto prefactor = s / bare_atomic_weight / bare_atomic_reweighting / double(nsamples);

    for (int i : range(nsamples)) {
      auto tau1 = data.tau_seg.get_random_pt(rng);
      auto tau2 = data.tau_seg.get_random_pt(rng);
      double dtau = double(tau2 - tau1);

      try {
        data.imp_trace.try_insert(tau1, op1_d);
        data.imp_trace.try_insert(tau2, op2_d);
        std::tie(atomic_weight, atomic_reweighting) = data.imp_trace.compute();
      } catch (rbt_insert_error const &) { atomic_weight = 0.; }

      data.imp_trace.cancel_insert();
      O_tau[closest_mesh_pt(dtau)] += prefactor * atomic_weight * atomic_reweighting;
    }
  }

  void measure_O_tau_ins::collect_results(mpi::communicator const &c) {
    O_tau        = mpi::all_reduce(O_tau, c);
    average_sign = mpi::all_reduce(average_sign, c);

    O_tau *= double(O_tau.mesh().size() - 1) / real(average_sign);

    // Assuming commuting operators so that O(0) = O(beta)
    // use that to counter the half-sized edge bind, by taking the average
    // on the boundary

    int last     = O_tau.mesh().size() - 1;
    auto average = O_tau[0] + O_tau[last];

    O_tau[0]    = average;
    O_tau[last] = average;
  }
} // namespace triqs_cthyb
