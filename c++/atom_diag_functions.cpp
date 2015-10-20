#include "./atom_diag_functions.hpp"

namespace cthyb {


double partition_function(atom_diag const& atom, double beta) {
 double z = 0;
 for (auto const& es : atom.get_eigensystem())
  for (auto e : es.eigenvalues) z += std::exp(-beta * e);
 return z;
}

// -----------------------------------------------------------------

block_gf<imtime> atomic_gf(atom_diag const& atom, double beta, std::map<std::string, indices_t> const& gf_struct,
                                      int n_tau, std::vector<std::pair<int, int>> const& excluded_states) {

 double z = partition_function(atom, beta);

 std::vector<std::string> block_names;
 std::vector<gf<imtime>> gf_blocks;
 auto const & fops = atom.get_fops();

 auto is_excluded = [&excluded_states](int A, int ia) {
  return std::find(excluded_states.begin(), excluded_states.end(), std::make_pair(A, ia)) != excluded_states.end();
 };

 for (auto const& block : gf_struct) {
  block_names.push_back(block.first);
  int bl_size = block.second.size();
  auto g = gf<imtime>{{beta, Fermion, n_tau}, {bl_size, bl_size}};

  for (int inner_index1 = 0; inner_index1 < bl_size; ++inner_index1)
   for (int inner_index2 = 0; inner_index2 < bl_size; ++inner_index2) {
    int n1 = fops[{block.first, block.second[inner_index1]}]; // linear_index of c
    int n2 = fops[{block.first, block.second[inner_index2]}]; // linear_index of c_dag

    for (int A = 0; A < atom.n_blocks(); ++A) {              // index of the A block. sum over all
     int B = atom.cdag_connection(n2, A);                // index of the block connected to A by operator c_n
     if (B == -1) continue;                             // no matrix element
     if (atom.c_connection(n1, B) != A) continue; //
     for (int ia = 0; ia < atom.get_block_dim(A); ++ia)
      for (int ib = 0; ib < atom.get_block_dim(B); ++ib) {
       auto Ea = atom.get_eigenvalue(A, ia);
       auto Eb = atom.get_eigenvalue(B, ib);
       if (is_excluded(A, ia) || is_excluded(B, ib)) continue;
       for (auto tau : g.mesh())
        g[tau](inner_index1, inner_index2) +=
            -atom.cdag_matrix(n2, A)(ib, ia) * atom.c_matrix(n1, B)(ia, ib) * std::exp(-(Eb - Ea) * tau - beta * Ea) / z;
      }
    }
   }
  g.singularity()(1) = 1.0;
  gf_blocks.push_back(std::move(g));
 }

 return make_block_gf(block_names, gf_blocks);
}

// -----------------------------------------------------------------

h_scalar_t trace_rho_op(block_matrix_t const& density_matrix, many_body_op_t const& op, atom_diag const& atom) {
 h_scalar_t result = 0;
 if (atom.n_blocks() != density_matrix.size()) TRIQS_RUNTIME_ERROR << "trace_rho_op : size mismatch : nmber of blocks differ";
 for (int bl = 0; bl < atom.n_blocks(); ++bl) {
  if (atom.get_block_dim(bl) != first_dim(density_matrix[bl]))
   TRIQS_RUNTIME_ERROR << "trace_rho_op : size mismatch : size of block " << bl << " differ";
  for (auto const& x : op) {
   auto b_m = atom.matrix_element_of_monomial(x.monomial, bl);
   if (b_m.first != -1) result += x.coef * dot_product(b_m.second, density_matrix[bl]);
  }
 }
 return result;
}

// -----------------------------------------------------------------

full_hilbert_space_state_t act(many_body_op_t const& op, full_hilbert_space_state_t const& st, atom_diag const& atom) {
 full_hilbert_space_state_t result(st.size());
 result() = 0;
 for (auto const& x : op) {
  for (int bl = 0; bl < atom.n_blocks(); ++bl) {
   auto b_m = atom.matrix_element_of_monomial(x.monomial, bl);
   // FIXME double : to be removed when ported to the right many_body_operator<double>
   if (b_m.first == -1) continue;
   result(atom.index_range_of_block(b_m.first)) += x.coef * b_m.second * st(atom.index_range_of_block(bl));
  }
 }
 return result;
}

//---------------------

std::vector<std::vector<double>> quantum_number_eigenvalues(many_body_op_t const& op, atom_diag const& atom) {

 auto commutator = op * atom.get_h_atomic() - atom.get_h_atomic() * op;
 if (!commutator.is_zero()) TRIQS_RUNTIME_ERROR << "The operator is not a quantum number";

 std::vector<std::vector<double>> result;

 for (int bl = 0; bl < atom.n_blocks(); ++bl) {
  auto dim = atom.get_block_dim(bl);
  result.push_back(std::vector<double>(dim, 0));
  for (auto const& x : op) {
   auto b_m = atom.matrix_element_of_monomial(x.monomial, bl);
   if (b_m.first != bl) continue;
   for (int i = 0; i < dim; ++i) result.back()[i] += x.coef * b_m.second(i, i);
  }
 }
 return result;
}

//---------------------

template<typename M> 
// require (ImmutableMatrix<M>)
bool is_diagonal(M const &m) {
 //auto r = trace(abs(m)); //0;
 //for (int i = 0; i < first_dim(m); ++i) r += abs(m(i, i));
 return ((sum(abs(m)) - trace(abs(m))) < 1.e-11);
}

std::vector<std::vector<double>> quantum_number_eigenvalues2(many_body_op_t const& op, atom_diag const& atom) {

 auto commutator = op * atom.get_h_atomic() - atom.get_h_atomic() * op;
 if (!commutator.is_zero()) TRIQS_RUNTIME_ERROR << "The operator is not a quantum number";

 auto d = atom.get_full_hilbert_space_dim();
 matrix<double> M(d, d);
 std::vector<std::vector<double>> result;

 for (int bl = 0; bl < atom.n_blocks(); ++bl) {
  auto dim = atom.get_block_dim(bl);
  for (auto const& x : op) {
   auto b_m = atom.matrix_element_of_monomial(x.monomial, bl);
   if (b_m.first == -1) continue;
   M(atom.index_range_of_block(b_m.first), atom.index_range_of_block(bl)) += x.coef * b_m.second;
   //for (int i = 0; i < dim; ++i) result.back()[i] += x.coef * b_m.second(i, i);
  }
 }
 // 
 if (!is_diagonal(M)) TRIQS_RUNTIME_ERROR << "The Matrix of the operator is not diagonal !!!";

 for (int bl = 0; bl < atom.n_blocks(); ++bl) {
   auto dim = atom.get_block_dim(bl);
   result.push_back(std::vector<double>(dim, 0));
   for (int i = 0; i < dim; ++i) result.back()[i] = M(atom.flatten_block_index(bl, i), atom.flatten_block_index(bl, i));
 }

 return result;
}
}
