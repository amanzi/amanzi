/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation of TreeVector, a nested, hierarchical data structure for PK
hiearchies.  This nested vector allows each physical PK to push back
Epetra_MultiVectors to store their solution, and allows MPCs to push back
TreeVectors in a nested format.  It also provides an implementation of the
Vector interface for use with time integrators/nonlinear solvers.
------------------------------------------------------------------------- */

#include "dbc.hh"
#include "TreeVector.hh"

namespace Amanzi {

TreeVector::TreeVector(std::string name) : name_(name) {};

TreeVector::TreeVector(std::string name, Teuchos::RCP<Epetra_MultiVector> &data) :
    name_(name) {
  data_.push_back(data);
};

TreeVector::TreeVector(std::string name, Teuchos::RCP<Epetra_Vector> &data) :
    name_(name) {
  data_.push_back(data);
};

TreeVector::TreeVector(std::string name, Teuchos::RCP<TreeVector> &subvec) : name_(name) {
  subvecs_.push_back(subvec);
};

TreeVector::TreeVector(std::string name,
                       std::vector< Teuchos::RCP<Epetra_MultiVector> > &dvec) :
    name_(name) {
  data_ = dvec;
};

TreeVector::TreeVector(std::string name, std::vector< Teuchos::RCP<TreeVector> > &subvecs) :
    name_(name) {
  subvecs_ = subvecs;
};

TreeVector& TreeVector::operator=(double value) {
  for (std::vector< Teuchos::RCP<Epetra_MultiVector> >::iterator datum = data_.begin();
       datum != data_.end(); ++datum) {
    (*datum)->PutScalar(value);
  }
  for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    *(*subvec) = value;
  }
  return *this;
};

TreeVector& TreeVector::operator=(const Epetra_Vector &value) {
  ASSERT(subvecs_.size() == 0);
  ASSERT(data_.size() == 1);
  Epetra_MultiVector datum = *data_[0];
  ASSERT(datum.NumVectors() == 1);
  *(datum(0)) = value;
  return *this;
};

TreeVector& TreeVector::operator=(const Epetra_MultiVector &value) {
  ASSERT(subvecs_.size() == 0);
  ASSERT(data_.size() == 1);

  *data_[0] = value;
  return *this;
};

void TreeVector::Scale(double value) {
  for (std::vector< Teuchos::RCP<Epetra_MultiVector> >::iterator datum = data_.begin();
       datum != data_.end(); ++datum) {
    (*datum)->Scale(value);
  }
  for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    (*subvec)->Scale(value);
  }
};

void TreeVector::Shift(double scalarA) {
  for (unsigned int i = 0; i != data_.size(); ++i) {
    Epetra_MultiVector work(*data_[i]);
    work.PutScalar(1.0);
    data_[i]->Update(scalarA, work, 1.0);
  }
  for (unsigned int i = 0; i != subvecs_.size(); ++i) {
    subvecs_[i]->Shift(scalarA);
  }
};

Vector& TreeVector::Update(double scalarA, const Vector& A, double scalarThis) {
  const TreeVector *A_ptr = dynamic_cast<const TreeVector *> (&A);
  ASSERT(A_ptr);
  ASSERT(A_ptr->subvecs_.size() == subvecs_.size());
  ASSERT(A_ptr->data_.size() == data_.size());

  for (unsigned int i = 0; i != data_.size(); ++i) {
    data_[i]->Update(scalarA, *(A_ptr->data_[i]), scalarThis);
  }
  for (unsigned int i = 0; i != subvecs_.size(); ++i) {
    subvecs_[i]->Update(scalarA, *A_ptr->subvecs_[i], scalarThis);
  }
  return *this;
};

Vector& TreeVector::Update(double scalarA, const Vector& A,
                   double scalarB, const Vector& B, double scalarThis) {
  const TreeVector *A_ptr = dynamic_cast<const TreeVector *> (&A);
  ASSERT(A_ptr);
  ASSERT(A_ptr->subvecs_.size() == subvecs_.size());
  ASSERT(A_ptr->data_.size() == data_.size());

  const TreeVector *B_ptr = dynamic_cast<const TreeVector *> (&B);
  ASSERT(B_ptr);
  ASSERT(B_ptr->subvecs_.size() == subvecs_.size());
  ASSERT(B_ptr->data_.size() == data_.size());

  for (unsigned int i = 0; i != data_.size(); ++i) {
    data_[i]->Update(scalarA, *(A_ptr->data_[i]), scalarB, *(B_ptr->data_[i]), scalarThis);
  }
  for (unsigned int i = 0; i != subvecs_.size(); ++i) {
    subvecs_[i]->Update(scalarA, *(A_ptr->subvecs_[i]), scalarB, *(B_ptr->subvecs_[i]),
                        scalarThis);
  }
  return *this;
};

int TreeVector::Dot(const Vector& other, double* result) const {
  const TreeVector *other_ptr = dynamic_cast<const TreeVector *> (&other);
  if (!other_ptr) return 1;
  if (other_ptr->subvecs_.size() != subvecs_.size()) return 2;
  if (other_ptr->data_.size() == data_.size()) return 3;

  for (unsigned int i = 0; i != data_.size(); ++i) {
    int ierr = data_[i]->Dot(*(other_ptr->data_[i]), result);
    if (ierr) return ierr;
  }
  for (unsigned int i = 0; i != subvecs_.size(); ++i) {
    int ierr = subvecs_[i]->Dot(*(other_ptr->subvecs_[i]), result);
    if (ierr) return ierr;
  }
  return 0;
};

int TreeVector::SubVector(std::string subname, Teuchos::RCP<TreeVector> &vec) {
  for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    if (subname == (*subvec)->Name()) {
      vec = *subvec;
      return 0;
    }
  }
  for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    int ierr = (*subvec)->SubVector(subname, vec);
    if (!ierr) return ierr;
  }
  return 1;
};

int TreeVector::SubVector(std::string subname,
                          Teuchos::RCP<const TreeVector> &vec) const {
  for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    if (subname == (*subvec)->Name()) {
      vec = *subvec;
      return 0;
    }
  }
  for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    int ierr = (*subvec)->SubVector(subname, vec);
    if (!ierr) return ierr;
  }
  return 1;
};

void TreeVector::PushBack(Teuchos::RCP<TreeVector>& subvec) {
  subvecs_.push_back(subvec);
};

void TreeVector::PushBack(Teuchos::RCP<Epetra_MultiVector>& data) {
  data_.push_back(data);
};

} // namespace
