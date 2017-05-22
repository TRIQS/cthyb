/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016, P. Seth, I. Krivenko, M. Ferrero and O. Parcollet
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

#include "./g4_l.hpp"

namespace cthyb {

  template <g4_channel Channel>
  measure_g4_l<Channel>::measure_g4_l(std::optional<g4_wll_t> &g4_wll_opt, qmc_data const &data, g4_measures_t &g4_measures)
     : data(data), g4_measures(g4_measures), average_sign(0) {

    const double beta = data.config.beta();

    order             = g4_measures.params.measure_g4_block_order;
    size_t n_l        = g4_measures.params.measure_g4_n_l;
    int n_bosonic     = g4_measures.params.measure_g4_n_bosonic;
    int nfft_buf_size = g4_measures.params.measure_g4_l_nfft_buf_size; 

    // Allocate the two-particle Green's function
    {
      gf_mesh<imfreq> mesh_w{beta, Boson, n_bosonic};
      gf_mesh<legendre> mesh_l{beta, Fermion, n_l};
      gf_mesh<cartesian_product<imfreq, legendre, legendre>> mesh_wll{mesh_w, mesh_l, mesh_l};

      g4_wll_opt = make_block2_gf(mesh_wll, g4_measures.gf_struct);
      g4_wll.rebind(*g4_wll_opt);
      g4_wll() = 0;
    }

    // Allocate the nfft buffers
    {
      nfft_buf.resize(mini_vector<int, 2>{g4_wll.size1(), g4_wll.size2()});

      gf_mesh<imfreq> mesh_w = std::get<0>(g4_wll(0, 0).mesh().components());

      for (auto const &m : g4_measures()) {
        auto s = m.target_shape;
        array<int, 6> buf_sizes(n_l, n_l, s[0], s[1], s[2], s[3]);
        buf_sizes() = nfft_buf_size;
        nfft_buf(m.b1.idx, m.b2.idx) = nfft_array_t<1, 6>{mesh_w, g4_wll(m.b1.idx, m.b2.idx).data(), buf_sizes};
      }
    }
  }

  template <g4_channel Channel> void measure_g4_l<Channel>::accumulate(mc_weight_t s) {

    s *= data.atomic_reweighting;
    average_sign += s;

    double beta = data.config.beta();
    int n_l     = std::get<1>(g4_wll(0, 0).mesh().components()).size();

    for (auto const &m : g4_measures()) {

      if (data.dets[m.b1.idx].size() == 0 || data.dets[m.b2.idx].size() == 0) continue;

      auto accumulate_impl = [&](op_t const &i, op_t const &j, op_t const &k, op_t const &l, det_scalar_t val) {

        tilde_p_gen p_l1_gen(beta), p_l2_gen(beta);
        double dtau = setup_times(p_l1_gen, p_l2_gen, i, j, k, l);

        for (int l1 : range(n_l)) {
          double p_l1 = p_l1_gen.next();
          for (int l2 : range(n_l)) {
            double p_l2 = p_l2_gen.next();
            mini_vector<int, 6> vec{l1, l2, i.second, j.second, k.second, l.second};
            nfft_buf(m.b1.idx, m.b2.idx).push_back({dtau}, vec, val * p_l1 * p_l2);
          }
        }
      };

      bool diag_block = (m.b1.idx == m.b2.idx);

      // Perform the accumulation looping over both determinants
      if (order == AABB || diag_block) {
        foreach (data.dets[m.b1.idx], [&](op_t const &i, op_t const &j, det_scalar_t M_ij) {
          foreach (data.dets[m.b2.idx], [&](op_t const &k, op_t const &l, det_scalar_t M_kl) {
            accumulate_impl(i, j, k, l, s * M_ij * M_kl); // Accumulate in legendre-nfft buffer
          })
            ;
        })
          ;
      }
      if (order == ABBA || diag_block) {
        foreach (data.dets[m.b1.idx], [&](op_t const &i, op_t const &l, det_scalar_t M_il) {
          foreach (data.dets[m.b2.idx], [&](op_t const &k, op_t const &j, det_scalar_t M_kj) {
            accumulate_impl(i, j, k, l, -s * M_il * M_kj); // Accumulate in legendre-nfft buffer
          })
            ;
        })
          ;
      }
    }
  }

  template <>
  double measure_g4_l<PH>::setup_times(tilde_p_gen &p_l1_gen, tilde_p_gen &p_l2_gen, op_t const &i, op_t const &j, op_t const &k, op_t const &l) {
    p_l1_gen.reset(i.first, j.first);
    p_l2_gen.reset(k.first, l.first);
    double dtau = 0.5 * double(i.first + j.first - k.first - l.first);
    return dtau;
  }

  template <>
  double measure_g4_l<PP>::setup_times(tilde_p_gen &p_l1_gen, tilde_p_gen &p_l2_gen, op_t const &i, op_t const &j, op_t const &k, op_t const &l) {
    p_l1_gen.reset(i.first, k.first);
    p_l2_gen.reset(j.first, l.first);
    double dtau = 0.5 * double(i.first + mult_by_int(j.first, 3) - mult_by_int(k.first, 3) - l.first);
    return dtau;
  }

  template <g4_channel Channel> void measure_g4_l<Channel>::collect_results(triqs::mpi::communicator const &c) {

    for (auto const &m : g4_measures()) { nfft_buf(m.b1.idx, m.b2.idx).flush(); }

    g4_wll = mpi_all_reduce(g4_wll, c);

    average_sign = mpi_all_reduce(average_sign, c);

    for (auto &g4_wll_block : g4_wll) {

      for (auto l : std::get<1>(g4_wll_block.mesh().components())) {
        auto _   = var_t{};
        double s = std::sqrt(2 * l + 1);
        g4_wll_block[_][l][_] *= s;
        g4_wll_block[_][_][l] *= s * (l % 2 ? 1 : -1);
      }

      g4_wll_block /= (real(average_sign) * data.config.beta());
    }
  }

  template class measure_g4_l<PP>;
  template class measure_g4_l<PH>;
}
