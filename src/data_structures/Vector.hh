/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#ifndef __VECTOR_HH__
#define __VECTOR_HH__

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

  // operations
  virtual void Scale(double value) = 0;
  //virtual void Shift(double shift) = 0;

  // this <- alpha*a + beta*b + gamma*this
  virtual Vector& Update(double scalarA, const Vector& A, double scalarThis) = 0;
  virtual Vector& Update(double scalarA, const Vector& A,
                         double scalarB, const Vector& B, double scalarThis) = 0;

  // returns the Vector (if this.name==subname) or sub-vector (if subname in
  // this) or throws exception
  virtual int SubVector(std::string subname, Teuchos::RCP<Vector>& subvec) = 0;
  virtual int SubVector(std::string subname, Teuchos::RCP<const Vector>& subvec) const = 0;



  // note that a Vector must also provide an innerProduct method, but it need
  // not work if the other Vector is not of the same type.  Implementations
  // likely use dynamic_cast?
  int Dot(const Vector& other, double* result);

};

} // namespace

#endif
