#incl
// FIXME : MERGE 
namespace triqs::operators {

  template <typename S> bool is_bosonic(many_body_operator_generic<S> const &op) {
    bool b = true, f = true;
    for (auto const &[monomial, coef] : op.get_monomials()) {
      b &= monomial.size() % 2 == 0;
      f &= monomial.size() % 2 == 1;
    }
    ALWAYS_EXPECTS(f != b, "Operator\n {} \n is neither bosonic nor fermionic", op);
    return b;
  }
} // namespace triqs::operators


