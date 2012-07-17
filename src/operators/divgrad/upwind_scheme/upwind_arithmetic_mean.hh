/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
// Author: Ethan Coon (ecoon@lanl.gov)
//
// Scheme for taking coefficients for div-grad operators from cells to
// faces.
// -----------------------------------------------------------------------------

#ifndef AMANZI_UPWINDING_ARITHMETICMEAN_SCHEME_
#define AMANZI_UPWINDING_ARITHMETICMEAN_SCHEME_

#include "upwinding.hh"

namespace Amanzi {
namespace Operators {

class UpwindArithmeticMean : public Upwinding {

public:

  UpwindArithmeticMean(std::string pkname,
                     std::string cell_coef,
                     std::string face_coef);

  void Update(const Teuchos::Ptr<State>& S);

  void CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const Teuchos::Ptr<CompositeVector>& face_coef);

private:

  std::string pkname_;
  std::string cell_coef_;
  std::string face_coef_;
};

} // namespace
} // namespace

#endif
