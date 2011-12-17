/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "dbc.hh"
#include "FlatVector.hh"

namespace Amanzi {

FlatVector::FlatVector(std::string name) : name_(name) {};

FlatVector::FlatVector(std::string name, Teuchos::RCP<Epetra_MultiVector> &vec) :
  name_(name),vec_(vec) {};

FlatVector::FlatVector(std::string name, Teuchos::RCP<Epetra_Vector> &vec) :
  name_(name),vec_(vec) {};

FlatVector& FlatVector::operator=(double value) {
  vec_->PutScalar(value);
  return *this;
};

FlatVector& FlatVector::operator=(const Epetra_Vector &value) {
  *vec_ = value;
};

FlatVector& FlatVector::operator=(const Epetra_MultiVector &value) {
  *vec_ = value;
};

void FlatVector::Scale(double value) {
  vec_->Scale(value);
};

FlatVector& FlatVector::Update(double scalarA, const Vector& A, double scalarThis) {
  const FlatVector *A_ptr = dynamic_cast<const FlatVector *> (&A);
  ASSERT(A_ptr);
  vec_->Update(scalarA, *(A_ptr->vec_), scalarThis);
  return *this;
};

FlatVector& FlatVector::Update(double scalarA, const Vector& A,
                        double scalarB, const Vector& B, double scalarThis) {
  const FlatVector *A_ptr = dynamic_cast<const FlatVector *> (&A);
  ASSERT(A_ptr);
  const FlatVector *B_ptr = dynamic_cast<const FlatVector *> (&B);
  ASSERT(B_ptr);

  vec_->Update(scalarA, *(A_ptr->vec_), scalarB, *(B_ptr->vec_), scalarThis);
  return *this;
};

int FlatVector::Dot(const Vector& other, double *result) const {
  const FlatVector *other_ptr = dynamic_cast<const FlatVector *> (&other);
  if (!other_ptr) return 1;

  return vec_->Dot(*(other_ptr->vec_), result);
};

} // namespace
