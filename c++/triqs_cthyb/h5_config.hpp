
friend void h5_write(h5::group g, op_desc const &op) {
  h5_write(g, "block", op.block_index);
  h5_write(g, "inner", op.inner_index);
  h5_write(g, "dagger", op.dagger);
}

// Writing of configuration out to a h5 for e.g. plotting
friend void h5_write(h5::group conf, std::string conf_group_name, configuration const &c) {
  h5::group conf_group = conf.create_group(conf_group_name);
  for (auto const &op : c) {
    // create group for given tau
    auto tau_group_name = std::to_string(double(op.first));
    h5::group tau_group = conf_group.create_group(tau_group_name);
    // in tau subgroup, write operator info
    h5_write(tau_group, op.second);
  }
}
