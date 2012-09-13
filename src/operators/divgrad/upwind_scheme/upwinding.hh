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

#ifndef AMANZI_UPWINDING_SCHEME_
#define AMANZI_UPWINDING_SCHEME_

namespace Amanzi {

// forward declaration
class State;

namespace Operators {

class Upwinding {

 public:
  virtual void
  Update(const Teuchos::Ptr<State>& S) = 0;
};

} // namespace
} // namespace

#endif
