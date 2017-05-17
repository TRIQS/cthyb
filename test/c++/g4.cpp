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

TEST(CtHyb, g4_measurments) {

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
  solver_core solver(beta, gf_struct, 1025, 2500);

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

  p.measure_g2_tau=true;
  p.measure_g2_n_tau=20;

  p.measure_g2_inu=true;  
  p.measure_g2_inu_fermionic=true;  
  p.measure_g2_n_inu=6;
    
  // Solve!
  solver.solve(p);

  std::cout << "--> solver done, now writing and reading the results.\n";
  
  // Save the results
  std::string filename = "g4";

  if (rank == 0) {
    triqs::h5::file G_file(filename + ".out.h5", 'w');
    h5_write(G_file, "G_tau", solver.G_tau());
    h5_write(G_file, "G2_tau", solver.G2_tau());
    h5_write(G_file, "G2_inu", solver.G2_inu());
  }

  if (rank == 0) {
    triqs::h5::file G_file(filename + ".ref.h5", 'r');

    g_tau_t g_tau;
    h5_read(G_file, "G_tau", g_tau);
    for( auto block_idx : range(g_tau.size()) )
      EXPECT_GF_NEAR(g_tau[block_idx], solver.G_tau()[block_idx]);

    g4_tau_t g4_tau;
    h5_read(G_file, "G2_tau", g4_tau);
    for( auto bidx1 : range(g4_tau.size1()) )
      for( auto bidx2 : range(g4_tau.size2()) )
	EXPECT_GF_NEAR(g4_tau(bidx1, bidx2), solver.G2_tau()(bidx1, bidx2));

    g4_iw_t g4_iw;
    h5_read(G_file, "G2_inu", g4_iw);
    for( auto bidx1 : range(g4_iw.size1()) )
      for( auto bidx2 : range(g4_iw.size2()) )
	EXPECT_GF_NEAR(g4_iw(bidx1, bidx2), solver.G2_inu()(bidx1, bidx2));
    
  }
}
MAKE_MAIN;
