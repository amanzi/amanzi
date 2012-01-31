/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Purely virtual interface for an ATS Vector.  Provides basic vector
functionality including scalar operations and an inner product.
Implementations of this interface may be used within ATS time integrators and
nonlinear solvers (eventually).
------------------------------------------------------------------------- */

#ifndef DATA_STRUCTURES_VECTOR_HH_
#define DATA_STRUCTURES_VECTOR_HH_

#include <string>

namespace Amanzi {

class Vector{
  // base class inherited by all Vectors used in Amanzi
public:

  // set data
  //  - guaranteed to not error
  virtual Vector& operator=(double value) = 0;

  //  - not guaranteed not to error
  virtual Vector& operator=(const Epetra_Vector &value) = 0;
  virtual Vector& operator=(const Epetra_MultiVector &value) = 0;

  // metadata
  virtual void SetName(std::string name) = 0;
  virtual std::string Name() = 0;

  // scalar operations
  virtual void Scale(double value) = 0;
  virtual void Shift(double value) = 0;

  // this <- alpha*a + beta*b + gamma*this
  virtual Vector& Update(double scalarA, const Vector& A, double scalarThis) = 0;
  virtual Vector& Update(double scalarA, const Vector& A,
                         double scalarB, const Vector& B, double scalarThis) = 0;

  virtual int Multiply (double scalarA, const Vector& A, const Vector& B, double scalarThis) =0;
  virtual int PutScalar(double scalar) = 0;
  virtual int NormInf(double* ninf) = 0;

  // note that a Vector must also provide an innerProduct method, but it need
  // not work if the other Vector is not of the same type.  Implementations
  // likely use dynamic_cast.
  virtual int Dot(const Vector& other, double* result) const = 0;
};
} // namespace

#endif
