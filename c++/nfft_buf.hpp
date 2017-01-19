/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017, N. Wentzell, I. Krivenko
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
#pragma once
#include <numeric>
#include <triqs/arrays.hpp>

#define NFFT_PRECISION_DOUBLE
#include <nfft3mp.h>

namespace cthyb {

using dcomplex = std::complex<double>;
using triqs::arrays::array_view;

// FIXME: remove when merged into TRIQS library master
template<typename T, int Rank>
T sum (triqs::arrays::mini_vector<T,Rank> const& a) {
 T res = {};
 for(int i = 0; i < Rank; ++i) res += a[i];
 return res;
}

template <int Rank> struct nfft_buf_t {

 // Possible future extensions:
 //  -Bosonic Matsubaras
 //  -Possibly remove tau shift and exponential in push_back

 /// Constructor
 nfft_buf_t(array_view<dcomplex, Rank> fiw_arr_, int buf_size_, double beta_, bool do_checks_ = false)
    : fiw_arr(fiw_arr_), buf_size(buf_size_), beta(beta_), do_checks(do_checks_) {

  // Capture frequency extents from fiw_arr and check that they are even (
  // i.e. fermionic matsubaras )
  auto freq_extents = triqs::arrays::get_shape(fiw_arr).to_vector(); // TODO Replace mini_vector with std::array
  std::vector<int> extents_int;
  for (int n : freq_extents) {
   if (n % 2 != 0)
    TRIQS_RUNTIME_ERROR << " dimension with uneven frequency count not "
                           "allowed in nfft_buf_t \n";
   extents_int.push_back(n);
  }

  // Init nfft_plan
  plan_ptr = std::make_unique<nfft_plan>();
  nfft_init(plan_ptr.get(), Rank, extents_int.data(), buf_size); // extents must be int here!
 }

 ~nfft_buf_t() {
  if (plan_ptr) nfft_finalize(plan_ptr.get());
 }

 /// Flush contents of the nfft buffer
 void flush() {
  // Trivial initialization of the remaining points
  for (int i = buf_counter; i < buf_size; ++i) {
   fx_arr()[i] = 0.0;
   for (int r = 0; r < Rank; ++r) x_arr()[i * Rank + r] = -0.5 + double(i) / buf_size;
  }
  do_nfft();
  buf_counter = 0;
 }

 // nfft_buffer needs to be uncopyable, because nfft_plan contains raw pointers
 nfft_buf_t(nfft_buf_t const &) = delete;
 nfft_buf_t(nfft_buf_t &&) = default;
 nfft_buf_t &operator=(nfft_buf_t const &) = delete;
 nfft_buf_t &operator=(nfft_buf_t &&) = default;

 /// Insert tau-vector {tau_1, tau_2, ... } \in [0,\beta)^Rank and
 /// corresponding f(tau) into the NFFT buffer
 void push_back(std::array<double, Rank> const &tau_arr, dcomplex ftau) {

  // Write the set of shifted and normalized tau values (i. e. x values) to
  // the NFFT buffer and sum taus
  double tau_sum = 0.0;
  for (int r = 0; r < Rank; ++r) {
   // Note: Nfft multi-arrays are stored in flattened arrays (c-order)
   x_arr()[buf_counter * Rank + r] = tau_arr[r] / beta - 0.5; // \in [-0.5, 0.5)
   tau_sum += tau_arr[r];                                     // Sum all tau values
  }

  // Write f(x) to nfft_plan-> The prefactor accounts for the Pi/beta offset
  // in fermionic Matsubaras
  fx_arr()[buf_counter] = std::exp(1_j * M_PI * tau_sum / beta) * ftau;

  ++buf_counter;

  // If buffer is full, perform transform
  if (is_full()) {
   do_nfft();
   buf_counter = 0;
  }
 }

 // Reset fiw_arr view
 void set_fiw_arr(array_view<dcomplex, Rank> new_fiw_arr) {
  using triqs::arrays::get_shape;
  assert(get_shape(new_fiw_arr) == get_shape(fiw_arr));
  fiw_arr.rebind(new_fiw_arr);
 }

 // Function to check whether buffer is empty
 bool is_empty() const { return buf_counter == 0; }

 // Function to check whether buffer is filled
 bool is_full() const { return (buf_counter >= buf_size); }

 private:
 // Triqs array to contain the final NFFT output in matsubara frequencies
 array_view<dcomplex, Rank> fiw_arr;

 // Number of tau points for the nfft
 int buf_size;

 // Inverse temperature of fiw_arr
 double beta;

 // Switch for testing in nfft
 bool do_checks;

 // Nfft3 plan that allocates memory and performs NFFT transform
 std::unique_ptr<nfft_plan> plan_ptr;

 // Counter for elements currently in buffer
 int buf_counter = 0;

 // Get pointer to array containing x values for the NFFT transform
 double *x_arr() { return plan_ptr->x; }

 // Get pointer to array containing f(x) values for the NFFT transform
 dcomplex *fx_arr() { return reinterpret_cast<dcomplex *>(plan_ptr->f); }

 // Get pointer to array containing the NFFT output h(k)
 const dcomplex *fk_arr() const { return reinterpret_cast<dcomplex *>(plan_ptr->f_hat); }

 // Perform NFFT transform and accumulate inside fiw_arr
 void do_nfft() {

  // NFFT Library precomputation and checks
  if (plan_ptr->flags & PRE_ONE_PSI) nfft_precompute_one_psi(plan_ptr.get());
  if (do_checks) { // Check validity of NFFT parameters
   const char *error_str = nfft_check(plan_ptr.get());
   if (error_str != 0) TRIQS_RUNTIME_ERROR << "Error in NFFT module: " << error_str << "\n";
  }

  // Execute transform
  nfft_adjoint(plan_ptr.get());

  // Accumulate results in fiw_arr. Care to normalize results afterwards
  int count = 0;
  for (auto fiw_itr = fiw_arr.begin(); fiw_itr != fiw_arr.end(); ++fiw_itr) {
   int factor = (sum(fiw_itr.indices()) % 2 ? -1 : 1);
   *fiw_itr += fk_arr()[count] * factor;
   ++count;
  }
 }
};
}
