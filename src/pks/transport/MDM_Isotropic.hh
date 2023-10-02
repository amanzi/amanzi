/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

  Isotropic mechanical dispersion model, primarily for code testing.
*/

#ifndef AMANZI_MDM_ISOTROPIC_HH_
#define AMANZI_MDM_ISOTROPIC_HH_

// TPLs
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "Factory.hh"
#include "Function.hh"
#include "Point.hh"
#include "Tensor.hh"

// Transport
#include "MDM.hh"

namespace Amanzi {
namespace Transport {

class MDM_Isotropic : public MDM {
 public:
  explicit MDM_Isotropic(Teuchos::ParameterList& plist);
  ~MDM_Isotropic(){};

  // Required methods from the base class
  // -- scalar dispersion tensor.
  WhetStone::Tensor
  mech_dispersion(double t,
                  const AmanziGeometry::Point& xc,
                  const AmanziGeometry::Point& u,
                  int axi_symmetry,
                  double wc,
                  double phi) const;


  // -- the model is valid if at least one parameter is not zero.
  bool is_valid() const { return (alpha_ != 0.0); }

 private:
  void ParseAlpha_(Teuchos::ParameterList& plist, const std::string& keyword);

 private:
  mutable double alpha_;
  std::unique_ptr<Function> alpha_func_;
  bool dispersivity_;

  static Utils::RegisteredFactory<MDM, MDM_Isotropic> reg_;
};

} // namespace Transport
} // namespace Amanzi

#endif
