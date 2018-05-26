#include <fstream>

inline bool exists(std::string filename) {
  std::ifstream f(filename.c_str());
  bool e = f.good();
  if (e) f.close();
  return e;
}
