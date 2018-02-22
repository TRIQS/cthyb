
#include "./G2_iw_acc.hpp"

namespace cthyb {
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
	  G2[w, n1, n2](i, j, k, l) << G2[w, n1, n2](i, j, k, l)
	    + s * M_ij[n1, n1 + w](i, j) * M_kl[n2 + w, n2](k, l);
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
	  G2[w, n1, n2](i, j, k, l) << G2[w, n1, n2](i, j, k, l)
          - s * M_il[n1, n2](i, l) * M_kj[n2 + w, n1 + w](k, j);
  }
  
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
} // namespace cthyb
