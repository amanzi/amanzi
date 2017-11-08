/* -*-  mode: c++; indent-tabs-mode: nil -*- */

// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
// Author: Ethan Coon (ecoon@lanl.gov)
//
// Scheme for taking coefficients for div-grad operators from cells to
// faces.
// -----------------------------------------------------------------------------

#ifndef AMANZI_UPWINDING_CELLCENTERED_SCHEME_
#define AMANZI_UPWINDING_CELLCENTERED_SCHEME_

#include "upwinding.hh"

namespace Amanzi {

class State;
class CompositeVector;

namespace Operators {

class UpwindCellCentered : public Upwinding {

public:

  UpwindCellCentered(std::string pkname,
                     std::string cell_coef,
                     std::string face_coef);

  void Update(const Teuchos::Ptr<State>& S,
              const Teuchos::Ptr<Debugger>& db=Teuchos::null);


  void CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const Teuchos::Ptr<CompositeVector>& face_coef);

  virtual std::string
  CoefficientLocation() { return "standard: cell"; }

private:

  std::string pkname_;
  std::string cell_coef_;
  std::string face_coef_;
};

} // namespace
} // namespace

#endif
