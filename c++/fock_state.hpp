/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by I. Krivenko, M. Ferrero, O. Parcollet
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

#ifndef TRIQS_CTQMC_KRYLOV_FOCK_STATE
#define TRIQS_CTQMC_KRYLOV_FOCK_STATE
#include <ostream>
#include <functional>

namespace triqs { namespace app { namespace impurity_solvers { namespace ctqmc_krylov {
//#define QUICK_APPLY_OP

#ifdef QUICK_APPLY_OP

typedef  uint64_t fock_state;

#else

 class fock_state {
  public:
      
      typedef uint64_t storage_t;
  
  private:
      
      union {
          storage_t all;
          uint8_t bytes[sizeof(storage_t)];
      } bitfield;
            
      static int count_table[256];
      
  protected:
      
      friend class complete_hilbert_space;
      
  //    fock_state(std::size_t i) : bitfield({i}) {}
  
  public:

   explicit fock_state(uint64_t i) : bitfield({i}) {}
   fock_state(fock_state const &) = default;
   fock_state(fock_state &&) = default;
   fock_state & operator = (fock_state const &) = default;

   operator storage_t() const { return bitfield.all; }

   bool filling_number(std::size_t mode) const { return bitfield.all & (static_cast<storage_t>(1) << mode); }

   int num_particles() const
   {
#ifndef OPTIMIZE_BITHACK
    int count = 0;
    for(int b=0; b<sizeof(storage_t); ++b) count += count_table[bitfield.bytes[b]];
    return count;
#else
    // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetNaive    
    // // option 1, for at most 14-bit values in v:
    return (bitfield.all * 0x200040008001ULL & 0x111111111111111ULL) % 0xf;
#endif
   }

   int num_particles_below(std::size_t mode) const
   {
    auto tmp_bitfield = bitfield ;
    tmp_bitfield.all &= (static_cast<storage_t>(1) << mode) - 1;
#ifndef OPTIMIZE_BITHACK
    int used_bytes = (mode+7)/8;
    int count = 0;
    for(int b=0; b<used_bytes; ++b) count += count_table[tmp_bitfield.bytes[b]];
    return count;
#else
    return (tmp_bitfield.all * 0x200040008001ULL & 0x111111111111111ULL) % 0xf;
#endif
   }

   void set_to_0(std::size_t mode) { bitfield.all &= ~(static_cast<storage_t>(1) << mode); }
   void set_to_1(std::size_t mode) { bitfield.all |= (static_cast<storage_t>(1) << mode); }

   friend std::ostream& operator<<(std::ostream& os, const fock_state& fs) { return os << fs.bitfield.all; }
 };

 int fock_state::count_table[]  =
 {
  0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
 };

#endif

}}}}

#ifndef QUICK_APPLY_OP

// hash function
namespace std {
 template<> class hash<triqs::app::impurity_solvers::ctqmc_krylov::fock_state> :
  public hash<triqs::app::impurity_solvers::ctqmc_krylov::fock_state::storage_t> {};
}
#endif


#endif
