/*
  This is the transport component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
  Usage: 
*/

#ifndef AMANZI_DIFFUSION_MODEL_HH_
#define AMANZI_DIFFUSION_MODEL_HH_

#include <vector>
#include <string>

#include "TransportDefs.hh"

namespace Amanzi {
namespace Transport {

class DiffusionModel {
 public:
  DiffusionModel() {};
  ~DiffusionModel() {};

  double FindComponentValue(const std::string tcc_name) {
    int nnames = names_.size();
    for (int i = 0; i < nnames; ++i) {
      if (names_[i] == tcc_name) return values_[i];
    }
    return 0.0;
  }

  // access
  std::vector<std::string>& names() { return names_; }
  std::vector<double>& values() { return values_; }

 private:
  std::vector<std::string> names_;
  std::vector<double> values_;
};

}  // namespace Transport
}  // namespace Amanzi

#endif

