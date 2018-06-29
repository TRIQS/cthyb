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

#include "./O_tau_ins.hpp"

namespace triqs_cthyb {

  using namespace triqs::gfs;

  measure_O_tau_ins::measure_O_tau_ins(std::optional<gf<imtime, scalar_valued>> &O_tau_opt, qmc_data const &data, int n_tau,
                                       many_body_op_t const &op1, many_body_op_t const &op2)
     : data(data), average_sign(0), op1(op1), op2(op2) {
    O_tau_opt = gf<imtime, scalar_valued>{{data.config.beta(), Fermion, n_tau}};
    O_tau.rebind(*O_tau_opt);
    O_tau() = 0.0;

    int aux_idx1 = data.imp_trace.add_aux_operator(op1);
    int aux_idx2 = data.imp_trace.add_aux_operator(op2);

    op1_d = op_desc{0, 0, true, aux_idx1};
    op2_d = op_desc{0, 0, true, aux_idx2};
  }

  void measure_O_tau_ins::accumulate(mc_weight_t s) {
    s *= data.atomic_reweighting;
    average_sign += s;

    double eps = 1e-14;
    double beta = O_tau.mesh().domain().beta;

    mc_weight_t trace_val;
    mc_weight_t bare_trace_val = data.imp_trace.compute().first;

    {

      double t1 = 0.;
      auto tau1 = data.tau_seg.make_time_pt(t1);

      for (auto t2 : O_tau.mesh()) {

        double eps = 0;
        if (t2 == 0.) eps = -1e-14;  // This should not be needed FIXME
        if (t2 == beta) eps = 1e-14; // This should not be needed FIXME

        auto tau2 = data.tau_seg.make_time_pt(t2 - eps);

        try {

          data.imp_trace.try_insert(tau1, op1_d);
          data.imp_trace.try_insert(tau2, op2_d);

          trace_val = data.imp_trace.compute().first;

        } catch (rbt_insert_error const &) {
          //std::cerr << "Insert error : recovering ... " << std::endl;
          trace_val = 0.;
        }

        data.imp_trace.cancel_insert();
        O_tau[t2] += trace_val / bare_trace_val;
      }
    }
  }

  void measure_O_tau_ins::collect_results(triqs::mpi::communicator const &c) {

    O_tau        = mpi_all_reduce(O_tau, c);
    average_sign = mpi_all_reduce(average_sign, c);
    O_tau /= -real(average_sign);
  }

} // namespace triqs_cthyb
