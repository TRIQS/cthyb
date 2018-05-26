#pragma once
#include <triqs/arrays.hpp>

namespace triqs {
  namespace arrays {

    //// FIXME TO BE MOVED INTO THE LIB AT MERGE

    /// <psi | M | psi> where M is by block and psi in a vector in the full Hilbert space.
    template <typename T> T vec_bdm_vec(vector<T> const &psi1, std::vector<matrix<T>> const &M, vector<T> const &psi2) {
      T r             = 0;
      int block_start = 0; // the index in the full Hilbert space of the first vector in the current block
      for (int bl = 0; bl < M.size(); ++bl) {
        int d = first_dim(M[bl]);
        if (d != second_dim(M[bl])) TRIQS_RUNTIME_ERROR << "vec_bdm_vec: the matrix of block " << bl << " is not square !";
        // not optimal: use psi1(range(block_start, block_start + a +1)) + BLAS if necessary FIXME
        for (int a = 0; a < d; ++a) {
          T tmp = 0;
          for (int b = 0; b < d; ++b) tmp += M[bl](a, b) * psi2(block_start + b);
          r += triqs::utility::conj(psi1(block_start + a)) * tmp;
        }
        block_start += d;
      }
      return r;
    }

    // Tr (A^* B) : Not optimal. Is blas possible ? FIXME
    template <typename T> T dot_product(matrix<T> const &a, matrix<T> const &b) {
      T r      = 0;
      int dim1 = first_dim(a), dim2 = second_dim(a);
      if ((dim1 != second_dim(b)) || (dim2 != first_dim(b))) TRIQS_RUNTIME_ERROR << "dot_product of matrices : size mismatch";
      for (int i = 0; i < dim1; ++i)
        for (int j = 0; j < dim2; ++j) r += triqs::utility::conj(a(i, j)) * b(j, i);
      return r;
    }
  }
}
