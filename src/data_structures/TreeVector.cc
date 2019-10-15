/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include "dbc.hh"
#include "TreeVector.hh"


namespace Amanzi {

TreeVector::TreeVector(const Teuchos::RCP<const TreeVectorSpace>& space,
                       InitMode mode)
{
  map_ = space;
  InitMap_(mode);
}

TreeVector::TreeVector(const TreeVector& other, Teuchos::DataAccess access,
                       InitMode mode)
{
  if (access == Teuchos::DataAccess::View) {
    Errors::Message message("TreeVector: View semantic not supported.");
    throw(message);
  }
  map_ = Teuchos::rcp(new TreeVectorSpace(*other.map_));
  InitMap_(mode);
  if (mode == InitMode::COPY) { *this = other; }
}


void
TreeVector::InitMap_(InitMode mode)
{
  if (mode != InitMode::NOALLOC && map_->Data() != Teuchos::null) {
    data_ = Teuchos::rcp(new CompositeVector(map_->Data()));
  }

  for (auto i : *map_) {
    subvecs_.emplace_back(Teuchos::rcp(new TreeVector(i, mode)));
  }
}


TreeVector&
TreeVector::operator=(const TreeVector& other)
{
  if (&other != this) {
    // Ensure the maps match.
    AMANZI_ASSERT(map_->SubsetOf(*other.map_));

    if (other.data_ != Teuchos::null) { *data_ = *other.data_; }

    for (unsigned int i = 0; i != subvecs_.size(); ++i) {
      *subvecs_[i] = *other.subvecs_[i];
    }
  }
  return *this;
};

void
TreeVector::putScalar(double scalar)
{
  // Set all data of this node and all child nodes to scalar.
  if (data_ != Teuchos::null) { data_->putScalar(scalar); }
  for (const auto& subvec : subvecs_) { subvec->putScalar(scalar); }
}

void
TreeVector::random()
{
  // Set all data of this node and all child nodes to random.
  if (data_ != Teuchos::null) { data_->random(); }
  for (const auto& subvec : subvecs_) { subvec->random(); }
}

double
TreeVector::normInf() const
{
  // Take the L_Inf norm of this.
  double ninf = 0.0;
  double ninf_loc;

  if (data_ != Teuchos::null) {
    ninf_loc = data_->normInf();
    if (ninf_loc > ninf) ninf = ninf_loc;
  }

  for (const auto& subvec : subvecs_) {
    ninf_loc = subvec->normInf();
    if (ninf_loc > ninf) ninf = ninf_loc;
  }
  return ninf;
};

double
TreeVector::norm1() const
{
  // Take the L_1 norm of this.
  double n1 = 0.0;
  double n1_loc;

  if (data_ != Teuchos::null) {
    n1_loc = data_->norm1();
    n1 += n1_loc;
  }

  for (const auto& subvec : subvecs_) {
    n1_loc = subvec->norm1();
    n1 += n1_loc;
  }
  return n1;
};

double
TreeVector::norm2() const
{
  // Take the L_2 norm of this.
  double n2 = 0.0;
  double n2_loc;

  if (data_ != Teuchos::null) {
    n2_loc = data_->norm2();
    n2 += pow(n2_loc, 2);
  }

  for (const auto& subvec : subvecs_) {
    n2_loc = subvec->norm2();
    n2 += pow(n2_loc, 2);
  }
  return sqrt(n2);
};

void
TreeVector::Print(std::ostream& os) const
{
  // Print data to ostream for this node and all children.
  if (data_ != Teuchos::null) data_->Print(os);

  for (std::vector<Teuchos::RCP<TreeVector>>::const_iterator subvec =
         subvecs_.begin();
       subvec != subvecs_.end();
       ++subvec) {
    (*subvec)->Print(os);
  }
};


// this <- abs(this)
void
TreeVector::abs(const TreeVector& other)
{
  if (data_ != Teuchos::null) { data_->abs(*other.data_); }
  int i = 0;
  for (const auto& subvec : subvecs_) {
    subvec->abs(*other.subvecs_[i]);
    i++;
  }
}


void
TreeVector::scale(double value)
{
  // this <- value*this
  if (data_ != Teuchos::null) { data_->scale(value); }

  for (const auto& subvec : subvecs_) { subvec->scale(value); }
};

// int TreeVector::Shift(double value) {
//   // this <- this + scalarA
//   int ierr = 0;
//   if (data_ != Teuchos::null) {
//     ierr = data_->Shift(value);
//     if (ierr) return ierr;
//   }

//   for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec =
//   subvecs_.begin();
//        subvec != subvecs_.end(); ++subvec) {
//     ierr = (*subvec)->Shift(value);
//     if (ierr) return ierr;
//   }
//   return ierr;
// };


// this <- element-wise reciprocal(this)
void
TreeVector::reciprocal(const TreeVector& other)
{
  // this <- value*this
  if (data_ != Teuchos::null) { data_->reciprocal(*other.data_); }

  int i = 0;
  for (const auto& subvec : subvecs_) {
    subvec->reciprocal(*other.subvecs_[i]);
    i++;
  }
}


double
TreeVector::dot(const TreeVector& other) const
{
  double result = 0.0;
  if (data_ != Teuchos::null) { result = data_->dot(*other.data_); }

  for (unsigned int i = 0; i != subvecs_.size(); ++i) {
    result += subvecs_[i]->dot(*other.subvecs_[i]);
  }
  return result;
};

// this <- scalarA*A + scalarThis*this
void
TreeVector::update(double scalarA, const TreeVector& A, double scalarThis)
{
  if (data_ != Teuchos::null) { data_->update(scalarA, *A.data_, scalarThis); }
  for (unsigned int i = 0; i != subvecs_.size(); ++i) {
    subvecs_[i]->update(scalarA, *A.subvecs_[i], scalarThis);
  }
};

void
TreeVector::update(double scalarA, const TreeVector& A, double scalarB,
                   const TreeVector& B, double scalarThis)
{
  if (data_ != Teuchos::null) {
    data_->update(scalarA, *A.data_, scalarB, *B.data_, scalarThis);
  }
  for (unsigned int i = 0; i != subvecs_.size(); ++i) {
    subvecs_[i]->update(
      scalarA, *A.subvecs_[i], scalarB, *B.subvecs_[i], scalarThis);
  }
};

void
TreeVector::elementWiseMultiply(double scalarAB, const TreeVector& A,
                                const TreeVector& B, double scalarThis)
{
  if (data_ != Teuchos::null) {
    data_->elementWiseMultiply(scalarAB, *A.data_, *B.data_, scalarThis);
  }

  for (unsigned int i = 0; i != subvecs_.size(); ++i) {
    subvecs_[i]->elementWiseMultiply(
      scalarAB, *A.subvecs_[i], *B.subvecs_[i], scalarThis);
  }
};

// int TreeVector::ReciprocalelementWiseMultiply(double scalarAB, const
// TreeVector& A,
//         const TreeVector& B, double scalarThis) {
//   //  AMANZI_ASSERT(map_->SubsetOf(*A.map_));
//   //  AMANZI_ASSERT(map_->SubsetOf(*B.map_));

//   int ierr = 0;
//   if (data_ != Teuchos::null) {
//     ierr = data_->ReciprocalelementWiseMultiply(scalarAB, *A.data_, *B.data_,
//     scalarThis); if (ierr) return ierr;
//   }

//   for (unsigned int i = 0; i != subvecs_.size(); ++i) {
//     ierr = subvecs_[i]->ReciprocalelementWiseMultiply(scalarAB,
//     *A.subvecs_[i], *B.subvecs_[i], scalarThis); if (ierr) return ierr;
//   }
//   return ierr;
// };

Teuchos::RCP<TreeVector>
TreeVector::SubVector(int index)
{
  // Get a pointer to the sub-vector by index
  if (index < subvecs_.size()) { return subvecs_[index]; }
  return Teuchos::null;
};

Teuchos::RCP<const TreeVector>
TreeVector::SubVector(int index) const
{
  // Get a pointer to the sub-vector by index
  if (index < subvecs_.size()) { return subvecs_[index]; }
  return Teuchos::null;
};

void
TreeVector::SetSubVector(int i, const Teuchos::RCP<TreeVector>& subvec)
{
  if (!subvec->getMap()->LocallySameAs(*getMap()->SubVector(i))) {
    Errors::Message message(
      "TreeVector: SetSubVector called with incompatible vector.");
    throw(message);
  }
  subvecs_[i] = subvec;
};

void
TreeVector::SetData(const Teuchos::RCP<CompositeVector>& data)
{
  if (!data->getMap()->LocallySameAs(*getMap()->Data())) {
    Errors::Message message(
      "TreeVector: SetData called with incompatible vector.");
    throw(message);
  }
  data_ = data;
};


int
TreeVector::getGlobalLength() const
{
  int total = 0;
  if (data_ != Teuchos::null) { total += data_->getGlobalLength(); }

  for (std::vector<Teuchos::RCP<TreeVector>>::const_iterator subvec =
         subvecs_.begin();
       subvec != subvecs_.end();
       ++subvec) {
    total += (*subvec)->getGlobalLength();
  }
  return total;
};

} // namespace Amanzi
