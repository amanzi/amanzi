/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

  Anisotropic mechanical dispersion model.
*/

#ifndef AMANZI_MDM_BEAR_HH_
#define AMANZI_MDM_BEAR_HH_

// TPLs
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "Factory.hh"
#include "Point.hh"
#include "Tensor.hh"

// Transport
#include "MDM.hh"

namespace Amanzi {
namespace Transport {

class MDM_Bear : public MDM {
 public:
  explicit MDM_Bear(Teuchos::ParameterList& plist);
  ~MDM_Bear(){};

  // Required methods from the base class
  // -- dispersion tensor of rank 2.
  WhetStone::Tensor
  mech_dispersion(const AmanziGeometry::Point& u, int axi_symmetry, double wc, double phi) const;

  // -- the model is valid if at least one parameter is not zero.
  bool is_valid() const { return (alphaL_ + alphaT_ != 0.0); }

 private:
  double alphaL_, alphaT_;

  static Utils::RegisteredFactory<MDM, MDM_Bear> factory_;
};

} // namespace Transport
} // namespace Amanzi

#endif
