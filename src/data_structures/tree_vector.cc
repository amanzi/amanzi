/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Implementation of TreeVector, a nested, hierarchical data structure for PK
   hierarchies.  This nested vector allows each physical PK to push back
   Epetra_MultiVectors to store their solution, and allows MPCs to push back
   TreeVectors in a nested format.  It also provides an implementation of the
   Vector interface for use with time integrators/nonlinear solvers.
   ------------------------------------------------------------------------- */

#include "dbc.hh"
#include "tree_vector.hh"

namespace Amanzi {

// Basic constructor of an empty TreeVector
TreeVector::TreeVector(std::string name) : name_(name) {};

TreeVector::TreeVector(std::string name,
                       const TreeVector& other,
                       ConstructMode mode) :
    name_(name) {
  if (other.data_ != Teuchos::null) {
    data_ = Teuchos::rcp(new CompositeVector(*other.data_, mode));
  }

  for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator other_subvec =
         other.subvecs_.begin(); other_subvec != other.subvecs_.end(); ++other_subvec) {
    Teuchos::RCP<Amanzi::TreeVector> new_subvec
      = Teuchos::rcp(new TreeVector((*other_subvec)->name_, **other_subvec, mode));
    subvecs_.push_back(new_subvec);
  }
};

TreeVector::TreeVector(const TreeVector& other, ConstructMode mode) :
    name_(std::string("")) {
  if (other.data_ != Teuchos::null) {
    data_ = Teuchos::rcp(new CompositeVector(*other.data_, mode));
  }

  for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator other_subvec =
         other.subvecs_.begin(); other_subvec != other.subvecs_.end(); ++other_subvec) {
    Teuchos::RCP<Amanzi::TreeVector> new_subvec
      = Teuchos::rcp(new TreeVector((*other_subvec)->name_, **other_subvec, mode));
    subvecs_.push_back(new_subvec);
  }
};

TreeVector& TreeVector::operator=(const TreeVector &other) {
  if (&other != this) {
    // Copy values at this node and all child nodes.
    ASSERT(subvecs_.size() == other.subvecs_.size());

    if (other.data_ != Teuchos::null) {
      *data_ = *other.data_;
    }

    for (unsigned int i = 0; i != subvecs_.size(); ++i) {
      *subvecs_[i] = *other.subvecs_[i];
    }
  }
  return *this;
};

int TreeVector::PutScalar(double scalar) {
  // Set all data of this node and all child nodes to scalar.
  int ierr = 0;
  if (data_ != Teuchos::null) {
    ierr = data_->PutScalar(scalar);
    if (ierr) return ierr;
  }
  for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    ierr = (*subvec)->PutScalar(scalar);
    if (ierr) return ierr;
  }
  return ierr;
}

int TreeVector::NormInf(double* ninf) const {
  // Take the L_Inf norm of this.
  if (ninf == NULL) return 1;
  if (data_ == Teuchos::null && subvecs_.size() == 0) return 1;

  int ierr = 0;
  *ninf = 0.0;
  double ninf_loc;

  if (data_ != Teuchos::null) {
    ierr = data_->NormInf(&ninf_loc);
    if (ierr) return ierr;
    if (ninf_loc > *ninf) *ninf = ninf_loc;
  }

  for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    ierr = (*subvec)->NormInf(&ninf_loc);
    if (ierr) return ierr;
    if (ninf_loc > *ninf) *ninf = ninf_loc;
  }
  return ierr;
};

int TreeVector::Norm1(double* n1) const {
  // Take the L_1 norm of this.
  if (n1 == NULL) return 1;
  if (data_ == Teuchos::null && subvecs_.size() == 0) return 1;

  int ierr = 0;
  *n1 = 0.0;
  double n1_loc;

  if (data_ != Teuchos::null) {
    ierr = data_->Norm1(&n1_loc);
    if (ierr) return ierr;
    *n1 += n1_loc;
  }

  for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    ierr = (*subvec)->Norm1(&n1_loc);
    if (ierr) return ierr;
    *n1 += n1_loc;
  }
  return ierr;
};

int TreeVector::Norm2(double* n2) const {
  // Take the L_2 norm of this.
  if (n2 == NULL) return 1;
  if (data_ == Teuchos::null && subvecs_.size() == 0) return 1;

  int ierr = 0;
  *n2 = 0.0;
  double n2_loc;

  if (data_ != Teuchos::null) {
    ierr = data_->Norm2(&n2_loc);
    if (ierr) return ierr;
    *n2 += pow(n2_loc,2);
  }

  for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    ierr = (*subvec)->Norm2(&n2_loc);
    if (ierr) return ierr;
    *n2 += pow(n2_loc,2);
  }
  *n2 = sqrt(*n2);
  return ierr;
};

void TreeVector::Print(ostream &os) const {
  // Print data to ostream for this node and all children.
  os << name_ << std::endl;

  if (data_ != Teuchos::null) data_->Print(os);

  for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    (*subvec)->Print(os);
  }
};

int TreeVector::Scale(double value) {
  // this <- value*this
  int ierr = 0;
  if (data_ != Teuchos::null) {
    ierr = data_->Scale(value);
    if (ierr) return ierr;
  }

  for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    ierr = (*subvec)->Scale(value);
    if (ierr) return ierr;
  }
  return ierr;
};

int TreeVector::Shift(double value) {
  // this <- this + scalarA
  int ierr = 0;
  if (data_ != Teuchos::null) {
    ierr = data_->Shift(value);
    if (ierr) return ierr;
  }

  for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    ierr = (*subvec)->Shift(value);
    if (ierr) return ierr;
  }
  return ierr;
};

int TreeVector::Dot(const TreeVector& other, double* result) const {
  // compute the dot product of all components of the tree vector
  // viewed as one flat vector
  if (result == NULL) return 1;
  if (data_ == Teuchos::null && subvecs_.size() == 0) return 1;
  ASSERT(subvecs_.size() == other.subvecs_.size());

  int ierr = 0;
  *result = 0.0;
  if (data_ != Teuchos::null) {
    ierr = data_->Dot(*other.data_, result);
    if (ierr) return ierr;
  }

  for (unsigned int i = 0; i != subvecs_.size(); ++i) {
    double intermediate_result;
    ierr = subvecs_[i]->Dot(*other.subvecs_[i], &intermediate_result);
    if (ierr) return ierr;
    *result += intermediate_result;
  }
  return ierr;
};

// this <- scalarA*A + scalarThis*this
TreeVector& TreeVector::Update(double scalarA, const TreeVector& A, double scalarThis) {
  ASSERT(subvecs_.size() == A.subvecs_.size());

  if (data_ != Teuchos::null) {
    data_->Update(scalarA, *A.data_, scalarThis);
  }
  for (unsigned int i = 0; i != subvecs_.size(); ++i) {
    subvecs_[i]->Update(scalarA, *A.subvecs_[i], scalarThis);
  }
  return *this;
};

TreeVector& TreeVector::Update(double scalarA, const TreeVector& A,
        double scalarB, const TreeVector& B, double scalarThis) {
  ASSERT(subvecs_.size() == A.subvecs_.size());
  ASSERT(subvecs_.size() == B.subvecs_.size());

  if (data_ != Teuchos::null) {
    data_->Update(scalarA, *A.data_, scalarB, *B.data_, scalarThis);
  }
  for (unsigned int i = 0; i != subvecs_.size(); ++i) {
    subvecs_[i]->Update(scalarA, *A.subvecs_[i], scalarB, *B.subvecs_[i], scalarThis);
  }
  return *this;
};

int TreeVector::Multiply(double scalarAB, const TreeVector& A, const TreeVector& B,
                         double scalarThis) {
  ASSERT(subvecs_.size() == A.subvecs_.size());
  ASSERT(subvecs_.size() == B.subvecs_.size());

  int ierr = 0;
  if (data_ != Teuchos::null) {
    ierr = data_->Multiply(scalarAB, *A.data_, *B.data_, scalarThis);
    if (ierr) return ierr;
  }

  for (unsigned int i = 0; i != subvecs_.size(); ++i) {
    ierr = subvecs_[i]->Multiply(scalarAB, *A.subvecs_[i], *B.subvecs_[i], scalarThis);
    if (ierr) return ierr;
  }
  return ierr;
};

Teuchos::RCP<TreeVector> TreeVector::SubVector(std::string subname) {
  // Get a pointer to the sub-vector "subname".
  for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    if (subname == (*subvec)->name()) {
      return *subvec;
    }
  }
  return Teuchos::null;
};

Teuchos::RCP<const TreeVector> TreeVector::SubVector(std::string subname) const {
  // Get a pointer to the sub-vector "subname".
  for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    if (subname == (*subvec)->name()) {
      return *subvec;
    }
  }
  return Teuchos::null;
};

void TreeVector::PushBack(Teuchos::RCP<TreeVector>& subvec) {
  subvecs_.push_back(subvec);
};

} // namespace
