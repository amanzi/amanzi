/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author Ethan Coon

Function from R^d to R^n.

------------------------------------------------------------------------- */

#ifndef AMANZI_VECTOR_FUNCTION_HH_
#define AMANZI_VECTOR_FUNCTION_HH_

namespace Amanzi {

class VectorFunction {
public:
  virtual VectorFunction* Clone() const = 0;

  virtual int size() const = 0;

  virtual double* operator() (const double* xt) const = 0;

};

} // namespace

#endif
