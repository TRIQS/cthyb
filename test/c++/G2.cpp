/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, H. U.R. Strand
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

#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/gfs.hpp>
#include <triqs/test_tools/gfs.hpp>

#include "solver_core.hpp"

using namespace cthyb;
using triqs::operators::c;
using triqs::operators::c_dag;
using triqs::operators::n;
using namespace triqs::gfs;
using indices_type = triqs::operators::indices_t;

TEST(CtHyb, G2_measurments) {

  std::cout << "Welcome to the CTHYB solver\n";

  // Initialize mpi
  int rank = triqs::mpi::communicator().rank();

  // Parameters
  double beta = 2.0;
  double U    = 0.0;
  double mu   = 2.0;

  double V1       = 2.0;
  double V2       = 5.0;
  double epsilon1 = 0.0;
  double epsilon2 = 4.0;

  // GF structure
  enum spin { up, down };
  std::map<std::string, indices_type> gf_struct{{"up", {0}}, {"down", {0}}};
  auto n_up   = n("up", 0);
  auto n_down = n("down", 0);

  // define operators
  auto H = U * n_up * n_down;

  // Construct CTQMC solver
  int n_iw  = 1025;
  int n_tau = 2500;
  int n_l   = 10;
  solver_core solver(beta, gf_struct, n_iw, n_tau, n_l);

  // Set G0
  triqs::clef::placeholder<0> om_;
  auto g0_iw = gf<imfreq>{{beta, Fermion}, {1, 1}};
  g0_iw(om_) << om_ + mu - V1 * V1 / (om_ - epsilon1) - V2 * V2 / (om_ - epsilon2);
  for (int bl = 0; bl < 2; ++bl) solver.G0_iw()[bl] = triqs::gfs::inverse(g0_iw);

  // Solve parameters
  int n_cycles      = 500;
  auto p            = solve_parameters_t(H, n_cycles);
  p.random_name     = "";
  p.random_seed     = 123 * rank + 567;
  p.max_time        = -1;
  p.length_cycle    = 100;
  p.n_warmup_cycles = 1000;
  p.move_double     = false;

  p.measure_G2_tau   = true;
  p.measure_G2_n_tau = 3;

  p.measure_G2_iw          = true;
  p.measure_G2_n_fermionic = 3;

  p.measure_G2_iw_ph     = true;
  p.measure_G2_iw_pp     = true;
  p.measure_G2_n_bosonic = 5;

  p.measure_G2_iwll_pp = true;
  p.measure_G2_n_l  = 3;

  // Solve!
  solver.solve(p);

  std::cout << "--> solver done, now writing and reading the results.\n";

  // Save the results
  std::string filename = "G2";

  if (rank == 0) {
    triqs::h5::file G_file(filename + ".out.h5", 'w');
    h5_write(G_file, "G2_tau", (*solver.G2_tau)(0, 1));
    if(solver.G2_iw) h5_write(G_file, "G2_iw", (*solver.G2_iw)(0, 1));
    if(solver.G2_iw_ph) h5_write(G_file, "G2_iw_ph", (*solver.G2_iw_ph)(0, 1));
    if(solver.G2_iw_ph) h5_write(G_file, "G2_iw_pp", (*solver.G2_iw_pp)(0, 1));
    if(solver.G2_iwll_pp) h5_write(G_file, "G2_iwll_pp", (*solver.G2_iwll_pp)(0, 1));
  }

  if (rank == 0) {
    triqs::h5::file G_file(filename + ".ref.h5", 'r');

    {
      G2_tau_t::g_t G2_tau;
      h5_read(G_file, "G2_tau", G2_tau);
      EXPECT_GF_NEAR(G2_tau, (*solver.G2_tau)(0, 1));
    }

    {
      G2_iw_t::g_t G2_iw;
      h5_read(G_file, "G2_iw", G2_iw);
      if(solver.G2_iw) EXPECT_GF_NEAR(G2_iw, (*solver.G2_iw)(0, 1));
    }

    {
      G2_iw_t::g_t G2_iw;
      h5_read(G_file, "G2_iw_ph", G2_iw);
      if(solver.G2_iw_ph) EXPECT_GF_NEAR(G2_iw, (*solver.G2_iw_ph)(0, 1));
    }

    {
      G2_iw_t::g_t G2_iw;
      h5_read(G_file, "G2_iw_pp", G2_iw);
      if(solver.G2_iw_pp) EXPECT_GF_NEAR(G2_iw, (*solver.G2_iw_pp)(0, 1));
    }

    {
      G2_iwll_t::g_t G2_iwll_pp;
      h5_read(G_file, "G2_iwll_pp", G2_iwll_pp);
      if(solver.G2_iwll_pp) EXPECT_GF_NEAR(G2_iwll_pp, (*solver.G2_iwll_pp)(0, 1));
    }
  }
}
MAKE_MAIN;
