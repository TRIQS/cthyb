
#include "./G2_iw_acc.hpp"

namespace triqs_cthyb {
  namespace G2_iw {

  namespace {

    // Index placeholders
    clef::placeholder<0> i;
    clef::placeholder<1> j;
    clef::placeholder<2> k;
    clef::placeholder<3> l;

    // Frequency placeholders
    clef::placeholder<4> w;
    clef::placeholder<5> n1;
    clef::placeholder<6> n2;
    clef::placeholder<7> n3;

  } // namespace

  // -- Particle-hole

  template <> void accumulate_impl_AABB<G2_channel::PH>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) {

    /*
    G2(w, n1, n2)
    (i, j, k, l) << G2(w, n1, n2)(i, j, k, l)                    //
          + s * M_ij(n1, n1 + w)(i, j) * M_kl(n2 + w, n2)(k, l); // sign in lhs in fft
    */

    auto const b_mesh = std::get<0>(G2.mesh());
    auto const f_mesh = std::get<1>(G2.mesh());

    for (auto const &w : b_mesh)
      for (auto const &n1 : f_mesh)
        for (auto const &n2 : f_mesh)
	  for (auto const i : range(G2.target_shape()[0]))
	  for (auto const j : range(G2.target_shape()[1]))
	  for (auto const k : range(G2.target_shape()[2]))
	  for (auto const l : range(G2.target_shape()[3]))
	  G2[w, n1, n2](i, j, k, l) += s * M_ij[n1, n1 + w](i, j) * M_kl[n2 + w, n2](k, l);

    //G2[w, n1, n2](i, j, k, l) << G2[w, n1, n2](i, j, k, l)
    //	    + s * M_ij[n1, n1 + w](i, j) * M_kl[n2 + w, n2](k, l);

  }

  template <> void accumulate_impl_ABBA<G2_channel::PH>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_il, M_t const &M_kj) {

    /*
    G2(w, n1, n2)
    (i, j, k, l) << G2(w, n1, n2)(i, j, k, l)                    //
          - s * M_il(n1, n2)(i, l) * M_kj(n2 + w, n1 + w)(k, j); // sign in lhs in fft
    */

    auto const b_mesh = std::get<0>(G2.mesh());
    auto const f_mesh = std::get<1>(G2.mesh());

    for (auto const &w : b_mesh)
      for (auto const &n1 : f_mesh)
        for (auto const &n2 : f_mesh)
	  for (auto const i : range(G2.target_shape()[0]))
	  for (auto const j : range(G2.target_shape()[1]))
	  for (auto const k : range(G2.target_shape()[2]))
	  for (auto const l : range(G2.target_shape()[3]))
	  G2[w, n1, n2](i, j, k, l) -= s * M_il[n1, n2](i, l) * M_kj[n2 + w, n1 + w](k, j);

    //G2[w, n1, n2](i, j, k, l) << G2[w, n1, n2](i, j, k, l)
    //     - s * M_il[n1, n2](i, l) * M_kj[n2 + w, n1 + w](k, j);
  }


// -------------------------------------------------------------------------------------
// W=0 and arbitrary size of target space
// -------------------------------------------------------------------------------------

    template <> void accumulate_impl_AABB_w0<G2_channel::PH>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) {

    // gaint speed by not doing matsubara frequency arithmetic ...

    /*

    auto const f_mesh = std::get<1>(G2.mesh());

    auto g2 = G2.data();
    auto m_ij = M_ij.data();
    auto m_kl = M_kl.data();

    for (auto const w1 : range(g2.shape()[1]))
    for (auto const w2 : range(g2.shape()[2]))
    for (auto const i : range(g2.shape()[3]))
    for (auto const j : range(g2.shape()[4]))
    for (auto const k : range(g2.shape()[5]))
    for (auto const l : range(g2.shape()[6]))
      g2(0, w1, w2, i, j, k, l) += s * m_ij(w1, w1, i, j) * m_kl(w2, w2, k, l);
    */

    auto const b_mesh = std::get<0>(G2.mesh());
    auto const f_mesh = std::get<1>(G2.mesh());

    for (auto const &w : b_mesh)
      for (auto const &n1 : f_mesh)
        for (auto const &n2 : f_mesh)
	  for (auto const i : range(G2.target_shape()[0]))
	  for (auto const j : range(G2.target_shape()[1]))
	  for (auto const k : range(G2.target_shape()[2]))
	  for (auto const l : range(G2.target_shape()[3]))
	  G2[w, n1, n2](i, j, k, l) += s * M_ij[n1, n1](i, j) * M_kl[n2, n2](k, l);
  }

  template <> void accumulate_impl_ABBA_w0<G2_channel::PH>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_il, M_t const &M_kj) {

    // gaint speed by not doing matsubara frequency arithmetic ...

    /*
    auto const f_mesh = std::get<1>(G2.mesh());

    auto g2 = G2.data();
    auto m_il = M_il.data();
    auto m_kj = M_kj.data();
    
    for (auto const w1 : range(g2.shape()[1]))
    for (auto const w2 : range(g2.shape()[2]))
    for (auto const i : range(g2.shape()[3]))
    for (auto const j : range(g2.shape()[4]))
    for (auto const k : range(g2.shape()[5]))
    for (auto const l : range(g2.shape()[6]))
      g2(0, w1, w2, i, j, k, l) -= s * m_il(w1, w2, i, l) * m_kj(w2, w1, k, j);
    */

    auto const b_mesh = std::get<0>(G2.mesh());
    auto const f_mesh = std::get<1>(G2.mesh());

    for (auto const &w : b_mesh)
      for (auto const &n1 : f_mesh)
        for (auto const &n2 : f_mesh)
	  for (auto const i : range(G2.target_shape()[0]))
	  for (auto const j : range(G2.target_shape()[1]))
	  for (auto const k : range(G2.target_shape()[2]))
	  for (auto const l : range(G2.target_shape()[3]))
	  G2[w, n1, n2](i, j, k, l) -= s * M_il[n1, n2](i, l) * M_kj[n2, n1](k, j);    
  }

    template <> void accumulate_impl_AABB_w0<G2_channel::PP>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) { TRIQS_RUNTIME_ERROR << "accumulate_impl_AABB_w0<G2_channel::PP> not implemented"; }
    template <> void accumulate_impl_ABBA_w0<G2_channel::PP>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) { TRIQS_RUNTIME_ERROR << "accumulate_impl_ABBA_w0<G2_channel::PP> not implemented";}
    
    template <> void accumulate_impl_AABB_w0<G2_channel::AllFermionic>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) { TRIQS_RUNTIME_ERROR << "accumulate_impl_AABB_w0<G2_channel::AllFermionic> not implemented"; }
    template <> void accumulate_impl_ABBA_w0<G2_channel::AllFermionic>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) { TRIQS_RUNTIME_ERROR << "accumulate_impl_ABBA_w0<G2_channel::AllFermionic> not implemented";}
    

// -------------------------------------------------------------------------------------
// W=0 and size 1 target space
// -------------------------------------------------------------------------------------
    
  template <> void accumulate_impl_AABB_opt_w0<G2_channel::PH>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) {

    // gaint speed by not doing matsubara frequency arithmetic ...

    /*
    auto const f_mesh = std::get<1>(G2.mesh());

    for (auto const &n1 : range(f_mesh.size()))
      for (auto const &n2 : range(f_mesh.size()))
	G2.data()(0, n1, n2, 0, 0, 0, 0) += s * M_ij.data()(n1, n1, 0, 0) * M_kl.data()(n2, n2, 0, 0);
    */
    auto const b_mesh = std::get<0>(G2.mesh());
    auto const f_mesh = std::get<1>(G2.mesh());

    for (auto const &w : b_mesh)
      for (auto const &n1 : f_mesh)
        for (auto const &n2 : f_mesh)
	  G2[w, n1, n2](0, 0, 0, 0) += s * M_ij[n1, n1](0, 0) * M_kl[n2, n2](0, 0);
    
  }

  template <> void accumulate_impl_ABBA_opt_w0<G2_channel::PH>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_il, M_t const &M_kj) {

    // gaint speed by not doing matsubara frequency arithmetic ...

    /*
    auto const f_mesh = std::get<1>(G2.mesh());

    for (auto const &n1 : range(f_mesh.size()))
      for (auto const &n2 : range(f_mesh.size()))
	G2.data()(0, n1, n2, 0, 0, 0, 0) -= s * M_il.data()(n1, n2, 0, 0) * M_kj.data()(n2, n1, 0, 0);
    */

    auto const b_mesh = std::get<0>(G2.mesh());
    auto const f_mesh = std::get<1>(G2.mesh());

    for (auto const &w : b_mesh)
      for (auto const &n1 : f_mesh)
        for (auto const &n2 : f_mesh)
	  G2[w, n1, n2](0, 0, 0, 0) -= s * M_il[n1, n2](0, 0) * M_kj[n2, n1](0, 0);
    
  }

    template <> void accumulate_impl_AABB_opt_w0<G2_channel::PP>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) { TRIQS_RUNTIME_ERROR << "accumulate_impl_AABB_opt_w0<G2_channel::PP> not implemented"; }
    template <> void accumulate_impl_ABBA_opt_w0<G2_channel::PP>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) { TRIQS_RUNTIME_ERROR << "accumulate_impl_ABBA_opt_w0<G2_channel::PP> not implemented";}
    
    template <> void accumulate_impl_AABB_opt_w0<G2_channel::AllFermionic>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) { TRIQS_RUNTIME_ERROR << "accumulate_impl_AABB_opt_w0<G2_channel::AllFermionic> not implemented"; }
    template <> void accumulate_impl_ABBA_opt_w0<G2_channel::AllFermionic>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) { TRIQS_RUNTIME_ERROR << "accumulate_impl_ABBA_opt_w0<G2_channel::AllFermionic> not implemented";}

// -------------------------------------------------------------------------------------
// General frequency but size 1 target space
// -------------------------------------------------------------------------------------

    
  template <> void accumulate_impl_AABB_opt<G2_channel::PH>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) {

    auto const b_mesh = std::get<0>(G2.mesh());
    auto const f_mesh = std::get<1>(G2.mesh());
    for (auto const &w : b_mesh)
      for (auto const &n1 : f_mesh)
	for (auto const &n2 : f_mesh)
	  G2[w, n1, n2](0, 0, 0, 0) += s * M_ij[n1, n1 + w](0, 0) * M_kl[n2 + w, n2](0, 0);
  }

  template <> void accumulate_impl_ABBA_opt<G2_channel::PH>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_il, M_t const &M_kj) {

    auto const b_mesh = std::get<0>(G2.mesh());
    auto const f_mesh = std::get<1>(G2.mesh());
    for (auto const &w : b_mesh)
      for (auto const &n1 : f_mesh)
	for (auto const &n2 : f_mesh)
	  G2[w, n1, n2](0, 0, 0, 0) += -s * M_il[n1, n2](0, 0) * M_kj[n2 + w, n1 + w](0, 0);   
  }
    
  // -- Particle-particle

  template <> void accumulate_impl_AABB<G2_channel::PP>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) {
    G2(w, n1, n2)
    (i, j, k, l) << G2(w, n1, n2)(i, j, k, l)                    //
          + s * M_ij(n1, w - n2)(i, j) * M_kl(w - n1, n2)(k, l); // sign in lhs in fft
  }

  template <> void accumulate_impl_ABBA<G2_channel::PP>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_il, M_t const &M_kj) {
    G2(w, n1, n2)
    (i, j, k, l) << G2(w, n1, n2)(i, j, k, l)                    //
          - s * M_il(n1, n2)(i, l) * M_kj(w - n1, w - n2)(k, j); // sign in lhs in fft
  }

  // -- Fermionic

  template <>
  void accumulate_impl_AABB<G2_channel::AllFermionic>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) {

    int size_ij = M_ij.target_shape()[0];
    int size_kl = M_kl.target_shape()[0];

    auto const iw_mesh = std::get<0>(G2.mesh());
    using mesh_point_t = typename decltype(iw_mesh)::mesh_point_t;

    for (auto const &n1 : iw_mesh)
      for (auto const &n2 : iw_mesh)
        for (auto const &n3 : iw_mesh) {
          mesh_point_t n4{iw_mesh, n1.index() + n3.index() - n2.index()};
          for (int i : range(size_ij))
            for (int j : range(size_ij))
              for (int k : range(size_kl))
                for (int l : range(size_kl)) G2[{n1, n2, n3}](i, j, k, l) += s * M_ij[{n2, n1}](j, i) * M_kl[{n4, n3}](l, k);
        }

    //G2(n1, n2, n3)(i, j, k, l) << G2(n1, n2, n3)(i, j, k, l) + s * M_ij(n2, n1)(j, i) * M_kl(n1 + n3 - n2, n3)(l, k);
  }

  template <>
  void accumulate_impl_ABBA<G2_channel::AllFermionic>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_il, M_t const &M_kj) {

    int size_il = M_il.target_shape()[0];
    int size_kj = M_kj.target_shape()[0];

    auto const iw_mesh = std::get<0>(G2.mesh());
    using mesh_point_t = typename decltype(iw_mesh)::mesh_point_t;
    
    for (auto const &n1 : iw_mesh)
      for (auto const &n2 : iw_mesh)
        for (auto const &n3 : iw_mesh) {
          mesh_point_t n4{iw_mesh, n1.index() + n3.index() - n2.index()};
          for (int i : range(size_il))
            for (int j : range(size_kj))
              for (int k : range(size_kj))
                for (int l : range(size_il)) G2[{n1, n2, n3}](i, j, k, l) += -s * M_il[{n4, n1}](l, i) * M_kj[{n2, n3}](j, k);
        }

    //G2(n1, n2, n3)(i, j, k, l) << G2(n1, n2, n3)(i, j, k, l) - s * M_il(n1 + n3 - n2, n1)(l, i) * M_kj(n2, n3)(j, k);
  }
    
    /*
    
    template <>
  inline void accumulate_impl_AABB<G2_channel::PH>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) {
    G2(w, n1, n2)
    (i, j, k, l) << G2(w, n1, n2)(i, j, k, l)                    //
          + s * M_ij(n1, n1 + w)(i, j) * M_kl(n2 + w, n2)(k, l); // sign in lhs in fft
  }

    template <>
  inline void accumulate_impl_ABBA<G2_channel::PH>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_il, M_t const &M_kj) {
    G2(w, n1, n2)
    (i, j, k, l) << G2(w, n1, n2)(i, j, k, l)                    //
          - s * M_il(n1, n2)(i, l) * M_kj(n2 + w, n1 + w)(k, j); // sign in lhs in fft
  }

  // -- Particle-particle

  template <>
  inline void accumulate_impl_AABB<G2_channel::PP>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij, M_t const &M_kl) {
    G2(w, n1, n2)
    (i, j, k, l) << G2(w, n1, n2)(i, j, k, l)                    //
          + s * M_ij(n1, w - n2)(i, j) * M_kl(w - n1, n2)(k, l); // sign in lhs in fft
  }

  template <>
  inline void accumulate_impl_ABBA<G2_channel::PP>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_il, M_t const &M_kj) {
    G2(w, n1, n2)
    (i, j, k, l) << G2(w, n1, n2)(i, j, k, l)                    //
          - s * M_il(n1, n2)(i, l) * M_kj(w - n1, w - n2)(k, j); // sign in lhs in fft
  }

  // -- Fermionic

  template <>
  inline void accumulate_impl_AABB<G2_channel::AllFermionic>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_ij,
                                                                            M_t const &M_kl) {

    int size_ij = M_ij.target_shape()[0];
    int size_kl = M_kl.target_shape()[0];

    auto const iw_mesh = std::get<0>(G2.mesh());
    using mesh_point_t = typename decltype(iw_mesh)::mesh_point_t;

    for (auto const &n1 : iw_mesh)
      for (auto const &n2 : iw_mesh)
        for (auto const &n3 : iw_mesh) {
          mesh_point_t n4{iw_mesh, n1.index() + n3.index() - n2.index()};
          for (int i : range(size_ij))
            for (int j : range(size_ij))
              for (int k : range(size_kl))
                for (int l : range(size_kl)) G2[{n1, n2, n3}](i, j, k, l) += s * M_ij[{n2, n1}](j, i) * M_kl[{n4, n3}](l, k);
        }

    //G2(n1, n2, n3)(i, j, k, l) << G2(n1, n2, n3)(i, j, k, l) + s * M_ij(n2, n1)(j, i) * M_kl(n1 + n3 - n2, n3)(l, k);
  }

  template <>
  inline void accumulate_impl_ABBA<G2_channel::AllFermionic>(G2_iw_t::g_t::view_type G2, mc_weight_t s, M_t const &M_il,
                                                                            M_t const &M_kj) {

    int size_il = M_il.target_shape()[0];
    int size_kj = M_kj.target_shape()[0];

    auto const iw_mesh = std::get<0>(G2.mesh());
    using mesh_point_t = typename decltype(iw_mesh)::mesh_point_t;

    for (auto const &n1 : iw_mesh)
      for (auto const &n2 : iw_mesh)
        for (auto const &n3 : iw_mesh) {
          mesh_point_t n4{iw_mesh, n1.index() + n3.index() - n2.index()};
          for (int i : range(size_il))
            for (int j : range(size_kj))
              for (int k : range(size_kj))
                for (int l : range(size_il)) G2[{n1, n2, n3}](i, j, k, l) += -s * M_il[{n4, n1}](l, i) * M_kj[{n2, n3}](j, k);
        }

    //G2(n1, n2, n3)(i, j, k, l) << G2(n1, n2, n3)(i, j, k, l) - s * M_il(n1 + n3 - n2, n1)(l, i) * M_kj(n2, n3)(j, k);
  }
    */


    
  } // namespace G2_iw
} // namespace triqs_cthyb
