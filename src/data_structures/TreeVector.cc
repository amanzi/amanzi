/*
  Data Structures

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Implementation of TreeVector, a nested, hierarchical data structure for PK
  hierarchies.  This nested vector allows each physical PK to push back
  Epetra_MultiVectors to store their solution, and allows MPCs to push back
  TreeVectors in a nested format.  It also provides an implementation of the
  Vector interface for use with time integrators/nonlinear solvers.
*/

#include "dbc.hh"
#include "TreeVector.hh"


namespace Amanzi {

// Basic constructors
TreeVector::TreeVector() :
    map_(Teuchos::rcp(new TreeVectorSpace()))
{}

TreeVector::TreeVector(const Comm_ptr_type& comm) :
    map_(Teuchos::rcp(new TreeVectorSpace(comm)))
{}

TreeVector::TreeVector(const TreeVectorSpace& space, InitMode mode)
{
  map_ = Teuchos::rcp(new TreeVectorSpace(space));
  InitMap_(mode);
  if (mode == INIT_MODE_ZERO) {
    PutScalar(0.);
  }
}

TreeVector::TreeVector(const Teuchos::RCP<TreeVectorSpace>& space, InitMode mode)
{
  map_ = space;
  InitMap_(mode);
  if (mode == INIT_MODE_ZERO) {
    PutScalar(0.);
  }
}

TreeVector::TreeVector(const TreeVector& other, InitMode mode)
{
  map_ = Teuchos::rcp(new TreeVectorSpace(*other.map_));
  InitMap_(mode);
  if (mode == INIT_MODE_ZERO) {
    PutScalar(0.);
  } else if (mode == INIT_MODE_COPY) {
    *this = other;
  }
}


void TreeVector::InitMap_(InitMode mode) {
  if (mode != INIT_MODE_NOALLOC &&
      map_->Data() != Teuchos::null) {
    data_ = Teuchos::rcp(new CompositeVector(*map_->Data()));
  }

  for (TreeVectorSpace::iterator i=map_->begin();
       i!=map_->end(); ++i) {
    InitPushBack_(Teuchos::rcp(new TreeVector(*i, mode)));
  }
}


TreeVector& TreeVector::operator=(const TreeVector &other) {
  if (&other != this) {
    // Ensure the maps match.
    AMANZI_ASSERT(map_->SubsetOf(*other.map_));

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

int TreeVector::PutScalarMasterAndGhosted(double scalar) {
  // Set all data of this node and all child nodes to scalar.
  int ierr = 0;
  if (data_ != Teuchos::null) {
    ierr = data_->PutScalarMasterAndGhosted(scalar);
    if (ierr) return ierr;
  }
  for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    ierr = (*subvec)->PutScalarMasterAndGhosted(scalar);
    if (ierr) return ierr;
  }
  return ierr;
}

int TreeVector::PutScalarGhosted(double scalar) {
  // Set all data of this node and all child nodes to scalar.
  int ierr = 0;
  if (data_ != Teuchos::null) {
    ierr = data_->PutScalarGhosted(scalar);
    if (ierr) return ierr;
  }
  for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    ierr = (*subvec)->PutScalarGhosted(scalar);
    if (ierr) return ierr;
  }
  return ierr;
}

int TreeVector::Random() {
  // Set all data of this node and all child nodes to random.
  int ierr = 0;
  if (data_ != Teuchos::null) {
    ierr = data_->Random();
    if (ierr) return ierr;
  }
  for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    ierr = (*subvec)->Random();
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

void TreeVector::Print(std::ostream &os, bool data_io) const {
  // Print data to ostream for this node and all children.
  if (data_ != Teuchos::null) data_->Print(os, data_io);

  for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    (*subvec)->Print(os, data_io);
  }
};


// this <- abs(this)
int TreeVector::Abs(const TreeVector& other) {
  // this <- value*this
  int ierr = 0;
  if (data_ != Teuchos::null) {
    ierr = data_->Abs(*other.data_);
    if (ierr) return ierr;
  }
  for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    ierr = (*subvec)->Abs(*other.subvecs_[subvec - subvecs_.begin()]);
    if (ierr) return ierr;
  }
  return ierr;
}


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


// this <- element-wise reciprocal(this)
int TreeVector::Reciprocal(const TreeVector& other) {
  // this <- value*this
  int ierr = 0;
  if (data_ != Teuchos::null) {
    ierr = data_->Reciprocal(*other.data_);
    if (ierr) return ierr;
  }
  for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    ierr = (*subvec)->Reciprocal(*other.subvecs_[subvec - subvecs_.begin()]);
    if (ierr) return ierr;
  }
  return ierr;
}


int TreeVector::Dot(const TreeVector& other, double* result) const {
  // compute the dot product of all components of the tree vector
  // viewed as one flat vector
  if (result == NULL) return 1;
  if (data_ == Teuchos::null && subvecs_.size() == 0) return 1;
  //  if (!map_->SameAs(*other.map_)) return 1;

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
  //  AMANZI_ASSERT(map_->SubsetOf(*A.map_));

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
  //  AMANZI_ASSERT(map_->SubsetOf(*A.map_));
  //  AMANZI_ASSERT(map_->SubsetOf(*B.map_));

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
  //  AMANZI_ASSERT(map_->SubsetOf(*A.map_));
  //  AMANZI_ASSERT(map_->SubsetOf(*B.map_));

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

int TreeVector::ReciprocalMultiply(double scalarAB, const TreeVector& A,
        const TreeVector& B, double scalarThis) {
  //  AMANZI_ASSERT(map_->SubsetOf(*A.map_));
  //  AMANZI_ASSERT(map_->SubsetOf(*B.map_));

  int ierr = 0;
  if (data_ != Teuchos::null) {
    ierr = data_->ReciprocalMultiply(scalarAB, *A.data_, *B.data_, scalarThis);
    if (ierr) return ierr;
  }

  for (unsigned int i = 0; i != subvecs_.size(); ++i) {
    ierr = subvecs_[i]->ReciprocalMultiply(scalarAB, *A.subvecs_[i], *B.subvecs_[i], scalarThis);
    if (ierr) return ierr;
  }
  return ierr;
};

Teuchos::RCP<TreeVector> TreeVector::SubVector(int index) {
  // Get a pointer to the sub-vector by index
  if (index < subvecs_.size()) {
    return subvecs_[index];
  }
  return Teuchos::null;
};

Teuchos::RCP<const TreeVector> TreeVector::SubVector(int index) const {
  // Get a pointer to the sub-vector by index
  if (index < subvecs_.size()) {
    return subvecs_[index];
  }
  return Teuchos::null;
};

void TreeVector::PushBack(const Teuchos::RCP<TreeVector>& subvec) {
  map_->PushBack(subvec->map_);
  InitPushBack_(subvec);
};

void TreeVector::InitPushBack_(const Teuchos::RCP<TreeVector>& subvec) {
  subvecs_.push_back(subvec);
};


void TreeVector::SetData(const Teuchos::RCP<CompositeVector>& data) {
  data_ = data;
  if (map_->Data() == Teuchos::null || !map_->Data()->SameAs(data->Map())) {
    map_->SetData(Teuchos::rcp(new CompositeVectorSpace(data->Map())));
  }
};


int TreeVector::GlobalLength() const {
  int total = 0;
  if (data_ != Teuchos::null) {
    total += data_->GlobalLength();
  }

  for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator subvec = subvecs_.begin();
       subvec != subvecs_.end(); ++subvec) {
    total += (*subvec)->GlobalLength();
  }
  return total;
};

} // namespace
