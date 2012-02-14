/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Purely virtual interface for an ATS Vector.  Provides basic vector
functionality including scalar operations and an inner product.
Implementations of this interface may be used within ATS time integrators and
nonlinear solvers (eventually through templating?).
------------------------------------------------------------------------- */

#ifndef DATA_STRUCTURES_VECTOR_HH_
#define DATA_STRUCTURES_VECTOR_HH_

#include <string>

namespace Amanzi {

class Vector{
  // base class inherited by all Vectors used in Amanzi
public:

  // set values
  virtual int PutScalar(double scalar) = 0;

  // norms
  virtual int NormInf(double* result) const;
  //virtual int Norm1(double* result) const;
  //virtual int Norm2(double* result) const;

  // scalar operations
  virtual int Scale(double value) = 0;
  virtual int Shift(double value) = 0;

  // The following operations must exist, but need not support other vectors
  // of different types.

  // virtual int Dot(const Vector& other, double* result) const = 0;
  // virtual Vector& Update(double scalarA, const Vector& A, double scalarThis) = 0;
  // virtual Vector& Update(double scalarA, const Vector& A,
  //                        double scalarB, const Vector& B, double scalarThis) = 0;
  // virtual int Multiply(double scalarA, const Vector& A, const Vector& B,
  //                      double scalarThis) = 0;

};
} // namespace

#endif
