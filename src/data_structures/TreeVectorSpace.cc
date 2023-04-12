/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   ATS

   Implementation of TreeVectorSpace, a factory for TreeVectors.
   ------------------------------------------------------------------------- */

//!
#include "errors.hh"
#include "CompositeSpace.hh"
#include "CompositeVector.hh"
#include "TreeVectorSpace.hh"

namespace Amanzi {

// Constructor from data
TreeVectorSpace::TreeVectorSpace(const Teuchos::RCP<const CompositeSpace>& cvfac)
  : data_(cvfac), comm_(cvfac->getComm())
{}

// Copy constructor
TreeVectorSpace::TreeVectorSpace(const TreeVectorSpace& other) : comm_(other.comm_)
{
  if (other.data_ != Teuchos::null) { data_ = Teuchos::rcp(new CompositeSpace(*other.data_)); }

  for (auto other_subvec : other) {
    Teuchos::RCP<TreeVectorSpace> new_subvec = Teuchos::rcp(new TreeVectorSpace(*other_subvec));
    subvecs_.push_back(new_subvec);
  }
};


// Checks equality
bool
TreeVectorSpace::isSameAs(const TreeVectorSpace& other) const
{
  if (data_ == Teuchos::null && other.data_ != Teuchos::null) return false;
  if (data_ != Teuchos::null && other.data_ == Teuchos::null) return false;
  if (data_ != Teuchos::null && other.data_ != Teuchos::null && !data_->isSameAs(*other.data_))
    return false;
  if (subvecs_.size() != other.subvecs_.size()) return false;
  for (int i = 0; i != subvecs_.size(); ++i)
    if (!subvecs_[i]->isSameAs(*other.subvecs_[i])) return false;
  return true;
}

bool
TreeVectorSpace::locallySameAs(const TreeVectorSpace& other) const
{
  if (data_ == Teuchos::null && other.data_ != Teuchos::null) return false;
  if (data_ != Teuchos::null && other.data_ == Teuchos::null) return false;
  if (data_ != Teuchos::null && other.data_ != Teuchos::null && !data_->locallySameAs(*other.data_))
    return false;
  if (subvecs_.size() != other.subvecs_.size()) return false;
  for (int i = 0; i != subvecs_.size(); ++i)
    if (!subvecs_[i]->locallySameAs(*other.subvecs_[i])) return false;
  return true;
}

// Checks subset
bool
TreeVectorSpace::isSubsetOf(const TreeVectorSpace& other) const
{
  if (data_ != Teuchos::null && other.data_ == Teuchos::null) return false;
  if (data_ != Teuchos::null && !data_->isSubsetOf(*other.data_)) return false;
  if (subvecs_.size() != other.subvecs_.size()) return false;
  for (int i = 0; i != subvecs_.size(); ++i)
    if (!subvecs_[i]->isSubsetOf(*other.subvecs_[i])) return false;
  return true;
}

Teuchos::RCP<const TreeVectorSpace>
TreeVectorSpace::getSubVector(int index) const
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

// debugging information
void
TreeVectorSpace::print(std::ostream& os) const
{
  int n(0);
  for (auto it = begin(); it != end(); ++it) {
    auto cvs = (*it)->getData();
    if (cvs == Teuchos::null) {
      os << "SubVector: " << n++ << std::endl;
      (*it)->print(os);
    } else {
      os << "SubVector (leaf): " << n++ << std::endl;
      for (auto kt = cvs->begin(); kt != cvs->end(); ++kt) {
        os << "  component: " << *kt << " (" << cvs->getNumVectors(*kt) << ")\n";
      }
    }
  }
};


Teuchos::RCP<TreeVectorSpace>
createTVSwithOneLeaf(const CompositeVectorSpace& cvs)
{
  auto tvs = Teuchos::rcp(new TreeVectorSpace());
  tvs->setData(cvs.CreateSpace());
  return tvs;
}


} // namespace Amanzi
