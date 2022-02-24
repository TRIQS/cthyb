
// a few details common to all moves

namespace triqs_cthyb {

  static const char *debug_config_print_start = "=============================================================";
  static const char *debug_config_print_end   = "";

  // add the histogram with name in the map histos, if histos.
  // Multiple insertion with the same name has no effect, cf std::map
  //
  inline histogram *add_histo(std::string const &name, histo_map_t *histos, double beta) {
    if (!histos) return nullptr;
    auto new_histo = histos->insert({name, {.0, beta, 100}});
    return &(new_histo.first->second);
  }


} // namespace triqs_cthyb
