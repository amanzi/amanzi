/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#ifndef AMANZI_FLOW_SPECIFIC_STORAGE_STANDARD_MODEL_HH_
#define AMANZI_FLOW_SPECIFIC_STORAGE_STANDARD_MODEL_HH_

#include <string>

#include "SpecificStorage.hh"

namespace Amanzi {
namespace Flow {

class SpecificStorage_Standard : public SpecificStorage {
 public:
  SpecificStorage_Standard(Teuchos::ParameterList& plist) : SpecificStorage(plist)
  {
    beta_f_ = plist.get<double>("fluid compressibility");
    beta_m_ = plist.get<double>("matrix compressibility");
    g_ = plist.get<double>("gravity");
    rho_ = plist.get<double>("fluid density");
  }

  virtual double Value(double porosity) const { return (porosity * beta_f_ + beta_m_) * g_ * rho_; }

 private:
  double beta_f_, beta_m_, g_, rho_;

  static Utils::RegisteredFactory<SpecificStorage, SpecificStorage_Standard> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
