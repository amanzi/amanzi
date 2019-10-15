/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include "errors.hh"
#include "CompositeSpace.hh"
#include "TreeVectorSpace.hh"

namespace Amanzi {

// Constructor from data
TreeVectorSpace::TreeVectorSpace(
  const Teuchos::RCP<const CompositeSpace>& cvfac)
  : data_(cvfac), comm_(cvfac->Comm())
{}

// Copy constructor
TreeVectorSpace::TreeVectorSpace(const TreeVectorSpace& other)
  : comm_(other.comm_)
{
  if (other.data_ != Teuchos::null) {
    data_ = Teuchos::rcp(new CompositeSpace(*other.data_));
  }

  for (auto other_subvec : other) {
    Teuchos::RCP<TreeVectorSpace> new_subvec =
      Teuchos::rcp(new TreeVectorSpace(*other_subvec));
    subvecs_.push_back(new_subvec);
  }
};


// Checks equality
bool
TreeVectorSpace::SameAs(const TreeVectorSpace& other) const
{
  if (data_ == Teuchos::null && other.data_ != Teuchos::null) return false;
  if (data_ != Teuchos::null && other.data_ == Teuchos::null) return false;
  if (data_ != Teuchos::null && other.data_ != Teuchos::null &&
      !data_->SameAs(*other.data_))
    return false;
  if (subvecs_.size() != other.subvecs_.size()) return false;
  for (int i = 0; i != subvecs_.size(); ++i)
    if (!subvecs_[i]->SameAs(*other.subvecs_[i])) return false;
  return true;
}

bool
TreeVectorSpace::LocallySameAs(const TreeVectorSpace& other) const
{
  if (data_ == Teuchos::null && other.data_ != Teuchos::null) return false;
  if (data_ != Teuchos::null && other.data_ == Teuchos::null) return false;
  if (data_ != Teuchos::null && other.data_ != Teuchos::null &&
      !data_->LocallySameAs(*other.data_))
    return false;
  if (subvecs_.size() != other.subvecs_.size()) return false;
  for (int i = 0; i != subvecs_.size(); ++i)
    if (!subvecs_[i]->LocallySameAs(*other.subvecs_[i])) return false;
  return true;
}

// Checks subset
bool
TreeVectorSpace::SubsetOf(const TreeVectorSpace& other) const
{
  if (data_ != Teuchos::null && other.data_ == Teuchos::null) return false;
  if (data_ != Teuchos::null && !data_->SubsetOf(*other.data_)) return false;
  if (subvecs_.size() != other.subvecs_.size()) return false;
  for (int i = 0; i != subvecs_.size(); ++i)
    if (!subvecs_[i]->SubsetOf(*other.subvecs_[i])) return false;
  return true;
}

// Teuchos::RCP<TreeVector>
// TreeVectorSpace::CreateVector() const {
//   Teuchos::RCP<TreeVector> result = Teuchos::rcp(new TreeVector());
//   if (data_ != Teuchos::null) {
//     result->set_data(data_->CreateVector());
//   }

//   for (std::vector< Teuchos::RCP<const TreeVectorSpace> >::const_iterator
//            sub=subvecs_.begin(); sub!=subvecs_.end(); ++sub) {
//     result->PushBack((*sub)->CreateVector());
//   }
//   return result;
// };


// Teuchos::RCP<TreeVector>
// TreeVectorSpace::CreateVector(bool ghosted) const {
//   Teuchos::RCP<TreeVector> result = Teuchos::rcp(new TreeVector());
//   if (data_ != Teuchos::null) {
//     result->set_data(data_->CreateVector(ghosted));
//   }

//   for (std::vector< Teuchos::RCP<const TreeVectorSpace> >::const_iterator
//            sub=subvecs_.begin(); sub!=subvecs_.end(); ++sub) {
//     result->PushBack((*sub)->CreateVector(ghosted));
//   }
//   return result;
// };


Teuchos::RCP<const TreeVectorSpace>
TreeVectorSpace::SubVector(int index) const
{
  // Get a pointer to the sub-vector by index
  if (index < subvecs_.size()) { return subvecs_[index]; }
  return Teuchos::null;
};

void
TreeVectorSpace::PushBack(const Teuchos::RCP<TreeVectorSpace>& subvec)
{
  subvecs_.push_back(subvec);
};


} // namespace Amanzi
