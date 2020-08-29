/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Isotropic mechanical dispersion model, primarily for code testing.
*/

#ifndef AMANZI_MDM_ISOTROPIC_HH_
#define AMANZI_MDM_ISOTROPIC_HH_

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

class MDM_Isotropic : public MDM {
 public:
  explicit MDM_Isotropic(Teuchos::ParameterList& plist);
  ~MDM_Isotropic() {};
  
  // Required methods from the base class
  // -- scalar dispersion tensor.
  WhetStone::Tensor mech_dispersion(
      const AmanziGeometry::Point& u, int axi_symmetric, double s, double phi) const;

  // -- the model is valid if at least one parameter is not zero.
  bool is_valid() const { return (alpha_ != 0.0); }

 private:
  double alpha_;
  bool dispersivity_;

  static Utils::RegisteredFactory<MDM, MDM_Isotropic> factory_;
};

}  // namespace Transport
}  // namespace Amanzi
 
#endif
