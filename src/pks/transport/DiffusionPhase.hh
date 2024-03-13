/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

*/

#ifndef AMANZI_DIFFUSION_PHASE_HH_
#define AMANZI_DIFFUSION_PHASE_HH_

#include <vector>
#include <string>

#include "TransportDefs.hh"

namespace Amanzi {
namespace Transport {

class DiffusionPhase {
 public:
  DiffusionPhase(){};
  ~DiffusionPhase(){};

  int FindDiffusionValue(const std::string& tcc_name, double* md)
  {
    for (int i = 0; i < names_.size(); ++i) {
      if (names_[i] == tcc_name) {
        *md = values_[i];
        return 0;
      }
    }

    *md = 0.0;
    return -1;
  }

  // access
  std::vector<std::string>& names() { return names_; }
  std::vector<double>& values() { return values_; }

 private:
  std::vector<std::string> names_;
  std::vector<double> values_;
};

} // namespace Transport
} // namespace Amanzi

#endif
