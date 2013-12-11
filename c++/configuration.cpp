#include "configuration.hpp"

namespace cthyb_krylov {

configuration::configuration(double beta_, sorted_spaces const& sosp, bool use_cutoff, double cutoff) : beta_(beta_) {

 // std::cout  << " Cutoff BS : "<< sosp.n_subspaces() << std::endl ;
 //std::cout  << " Cutoff BS : "<< sosp << std::endl ;
 // Atomic partition function
 double z = 0;
 for (auto const& es : sosp.get_eigensystems())
  for (auto e : es.eigenvalues) {
   z += exp(-beta_ * e);
  }

 // Fill boundary_block_states with states which have significant weights
 int total_states = 0;
 for (int nsp = 0; nsp < sosp.n_subspaces(); ++nsp) {
  auto const& eigensystem = sosp.get_eigensystems()[nsp];
  for (size_t n = 0; n < eigensystem.eigenvalues.size(); ++n) {
   double prob = exp(-beta_ * eigensystem.eigenvalues[n]) / z;
   if (prob > cutoff) {
    boundary_block_states_ids.push_back(std::make_pair(nsp, n));
    total_states++;
   }
  }
 }

 // not true any more
 //std::cout << "Summing over " << total_states << " states in the outer trace ";
 //std::cout << "(the probability cutoff is set to " << cutoff << ")" << std::endl;
}

//-------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& out, configuration const& c) {
 out << "boundary_block_states_ids:" << std::endl;
 for (auto const& st : c.boundary_block_states_ids) out << "(" << st.first << "," << st.second << ") ";
 out << std::endl;
 for (auto const& op : c.oplist)
  out << "tau = " << op.first << " : " << (op.second.dagger ? "Cdag(" : "C(") << op.second.block_index << ","
      << op.second.inner_index << ")\n";
 return out;
}
}
