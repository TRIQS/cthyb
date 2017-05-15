#pragma once

#include <triqs/gfs.hpp>
#include <triqs/utility/variant_int_string.hpp>

namespace std {
  inline string to_string(triqs::utility::variant_int_string const &vis) { return triqs::utility::to_string(vis); }
}

namespace triqs {
  namespace gfs {

    /// The structure of the gf : block_name -> [...]= list of indices (int/string). FIXME Change to pair of vec<str> and vec<int> or vec<pair<str,int>>
    using block_gf_structure_t = std::map<std::string, std::vector<triqs::utility::variant_int_string>>;

    // Function template for block_gf initialization
    template <typename Val_t, typename Var_t> block_gf<Var_t, Val_t> make_block_gf(gf_mesh<Var_t> const &m, block_gf_structure_t const &gf_struct) {

      std::vector<gf<Var_t, Val_t>> gf_vec;
      std::vector<std::string> block_names;

      //for (auto const & [ bname, idx_lst ] : gf_struct) { // C++17
      for (auto const &bl : gf_struct) {
        auto &bname  = bl.first;
        auto bl_size = bl.second.size();
        block_names.push_back(bname);
        std::vector<std::string> indices;
        for (auto const &var : bl.second) apply_visitor([&indices](auto &&arg) { indices.push_back(std::to_string(arg)); }, var);
        gf_vec.emplace_back(m, make_shape(bl_size, bl_size), std::vector<std::vector<std::string>>{indices, indices});
      }

      return make_block_gf(std::move(block_names), std::move(gf_vec));
    }

    // default to matrix_valued gf
    template <typename Var_t> block_gf<Var_t, matrix_valued> make_block_gf(gf_mesh<Var_t> const &m, block_gf_structure_t const &gf_struct) {
      return make_block_gf<matrix_valued, Var_t>(m, gf_struct);
    }

    template <typename Var_t> block2_gf<Var_t, tensor_valued<4>> make_block2_gf(gf_mesh<Var_t> const &m, block_gf_structure_t const &gf_struct) {

      std::vector<std::vector<gf<Var_t, tensor_valued<4>>>> gf_vecvec;
      std::vector<std::string> block_names;

      for (auto const &bl1 : gf_struct) {
        auto &bname  = bl1.first;
        int bl1_size = bl1.second.size();
        block_names.push_back(bname);
        std::vector<std::string> indices1;
        for (auto const &var : bl1.second) apply_visitor([&indices1](auto &&arg) { indices1.push_back(std::to_string(arg)); }, var);

        std::vector<gf<Var_t, tensor_valued<4>>> gf_vec;
        for (auto const &bl2 : gf_struct) {
          int bl2_size = bl2.second.size();
          std::vector<std::string> indices2;
          for (auto const &var : bl2.second) apply_visitor([&indices2](auto &&arg) { indices2.push_back(std::to_string(arg)); }, var);
          gf_vec.emplace_back(m, make_shape(bl1_size, bl1_size, bl2_size, bl2_size),
                              std::vector<std::vector<std::string>>{indices1, indices1, indices2, indices2});
        }
        gf_vecvec.emplace_back(std::move(gf_vec));
      }
      return make_block2_gf(block_names, block_names, std::move(gf_vecvec));
    }

  } // namespace gfs
} // namespace triqs
