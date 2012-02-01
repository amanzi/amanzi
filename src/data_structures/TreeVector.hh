/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for TreeVector, a nested, hierarchical data structure for PK
hiearchies.  This nested vector allows each physical PK to push back
Epetra_MultiVectors to store their solution, and allows MPCs to push back
TreeVectors in a nested format.  It also provides an implementation of the
Vector interface for use with time integrators/nonlinear solvers.
------------------------------------------------------------------------- */

#ifndef DATA_STRUCTURES_TREEVECTOR_HH_
#define DATA_STRUCTURES_TREEVECTOR_HH_

#include <string>
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Vector.hh"

namespace Amanzi {

class TreeVector : public Vector {

public:
  // constructor the other three
  explicit TreeVector(std::string name);
  TreeVector(std::string name, const Teuchos::RCP<Epetra_MultiVector>&);
  TreeVector(std::string name, const Teuchos::RCP<Epetra_Vector>&);
  TreeVector(std::string name, const Teuchos::RCP<TreeVector>&);
  TreeVector(std::string name, const std::vector< Teuchos::RCP<Epetra_MultiVector> >&);
  TreeVector(std::string name, const std::vector< Teuchos::RCP<TreeVector> >&);
  TreeVector(std::string name, const TreeVector&);
  TreeVector(const TreeVector&);
  TreeVector(const Teuchos::RCP<TreeVector>&);

  // set data
  //  - guaranteed to not error
  TreeVector& operator=(double value);

  //  - not guaranteed not to error
  TreeVector& operator=(const Epetra_Vector &value);
  TreeVector& operator=(const Epetra_MultiVector &value);
  TreeVector& operator=(const TreeVector &value);
  

  // metadata
  void SetName(std::string name) { name_ = name; }
  std::string Name() { return name_; }

  // operations
  void Scale(double value);
  void Shift(double scalarA);
  int Dot(const Vector& other, double* result) const;

  // this <- scalarA*A + scalarThis*this
  Vector& Update(double scalarA, const Vector& A, double scalarThis);

  // this <- scalarA*A + scalarB*B + scalarThis*this
  Vector& Update(double scalarA, const Vector& A,
                 double scalarB, const Vector& B, double scalarThis);
  
  // this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
  int Multiply(double scalarAB, const Vector& A, const Vector& B, double scalarThis);
  
  // this <- scalar
  int PutScalar(double scalar);
  
  // return || this ||_{inf}
  int NormInf(double* ninf);

  // non-inherited extras

  void Print(ostream &os) const;

  int SubVector(std::string subname, Teuchos::RCP<TreeVector>& subvec);
  int SubVector(std::string subname, Teuchos::RCP<const TreeVector>& subvec) const;

  Teuchos::RCP<Epetra_MultiVector> operator[](int vecnum) { return data_[vecnum]; }
  Teuchos::RCP<const Epetra_MultiVector> operator[](int vecnum) const { return data_[vecnum]; }

  void PushBack(Teuchos::RCP<TreeVector>& subvec);
  void PushBack(Teuchos::RCP<Epetra_MultiVector>& data);

private:
  std::string name_;
  std::vector< Teuchos::RCP<Epetra_MultiVector> > data_;
  std::vector< Teuchos::RCP<TreeVector> > subvecs_;
};

} // namespace


#endif
