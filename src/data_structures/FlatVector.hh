/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#ifndef __FLATVECTOR_HH__
#define __FLATVECTOR_HH__

#include <string>
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Vector.hh"

namespace Amanzi {

class FlatVector : public Vector {

public:
  // constructor the other three
  explicit FlatVector(std::string name);
  FlatVector(std::string name, Teuchos::RCP<Epetra_MultiVector>&);
  FlatVector(std::string name, Teuchos::RCP<Epetra_Vector>&);

  // set data
  //  - guaranteed to not error
  FlatVector& operator=(double value);

  //  - not guaranteed not to error
  FlatVector& operator=(const Epetra_Vector &value);
  FlatVector& operator=(const Epetra_MultiVector &value);

  // metadata
  void SetName(std::string name) { name_ = name; }
  std::string Name() { return name_; }

  // operations
  void Scale(double value);
  //void shift(double shift);

  // this <- alpha*a + beta*b + gamma*this
  FlatVector& Update(double scalarA, const Vector& A, double scalarThis);
  FlatVector& Update(double scalarA, const Vector& A,
              double scalarB, const Vector& B, double scalarThis);

  // returns the Vector (if this.name==subname) or sub-vector (if subname in
  // this)
  int SubVector(std::string subname, Teuchos::RCP<Vector> subvec) { return 1; }

  int Dot(const Vector& other, double* result) const;

private:
  std::string name_;
  Teuchos::RCP<Epetra_MultiVector> vec_;
};

} // namespace


#endif
