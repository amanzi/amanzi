/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#ifndef __TREEVECTOR_HH__
#define __TREEVECTOR_HH__

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
  TreeVector(std::string name, Teuchos::RCP<Epetra_MultiVector>&);
  TreeVector(std::string name, Teuchos::RCP<Epetra_Vector>&);
  TreeVector(std::string name, Teuchos::RCP<Vector>&);
  TreeVector(std::string name, std::vector< Teuchos::RCP<Epetra_MultiVector> >&);
  TreeVector(std::string name, std::vector< Teuchos::RCP<Vector> >&);

  // set data
  //  - guaranteed to not error
  TreeVector& operator=(double value);

  //  - not guaranteed not to error
  TreeVector& operator=(const Epetra_Vector &value);
  TreeVector& operator=(const Epetra_MultiVector &value);

  // metadata
  void SetName(std::string name) { name_ = name; }
  std::string Name() { return name_; }

  // operations
  void Scale(double value);
  //void shift(double shift);

  // this <- alpha*a + beta*b + gamma*this
  TreeVector& Update(double scalarA, const Vector& A, double scalarThis);
  TreeVector& Update(double scalarA, const Vector& A,
                     double scalarB, const Vector& B, double scalarThis);

  // returns the Vector (if this.name==subname) or sub-vector (if subname in
  // this)
  int SubVector(std::string subname, Teuchos::RCP<Vector>& subvec);
  int SubVector(std::string subname, Teuchos::RCP<const Vector>& subvec) const;

  int Dot(const Vector& other, double* result) const;

private:
  std::string name_;
  std::vector< Teuchos::RCP<Epetra_MultiVector> > data_;
  std::vector< Teuchos::RCP<Vector> > subvecs_;
};

} // namespace


#endif
