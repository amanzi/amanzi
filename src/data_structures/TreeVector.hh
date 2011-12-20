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
  TreeVector(std::string name, Teuchos::RCP<TreeVector>&);
  TreeVector(std::string name, std::vector< Teuchos::RCP<Epetra_MultiVector> >&);
  TreeVector(std::string name, std::vector< Teuchos::RCP<TreeVector> >&);
  TreeVector(const TreeVector&);

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
  void Shift(double scalarA);
  int Dot(const Vector& other, double* result) const;

  // this <- alpha*a + beta*b + gamma*this
  Vector& Update(double scalarA, const Vector& A, double scalarThis);
  Vector& Update(double scalarA, const Vector& A,
                 double scalarB, const Vector& B, double scalarThis);

  // non-inherited extras
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
