/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon (ecoon@lanl.gov)

Wraps a CompositeMatrix as a TreeMatrix.
------------------------------------------------------------------------- */

#ifndef AMANZI_COMPOSITEMATRIX_AS_TREEMATRIX_HH_
#  define AMANZI_COMPOSITEMATRIX_AS_TREEMATRIX_HH_

#  include "Teuchos_RCP.hpp"

#  include "CompositeVector.hh"
#  include "CompositeVectorSpace.hh"
#  include "composite_matrix.hh"
#  include "TreeVector.hh"
#  include "TreeVectorSpace.hh"
#  include "tree_matrix.hh"

namespace Amanzi {

class CompositeMatrixAsTreeMatrix : public TreeMatrix {
 public:
  CompositeMatrixAsTreeMatrix(const Teuchos::RCP<const CompositeMatrix>& cm) : cm_(cm) {}

  CompositeMatrixAsTreeMatrix(const CompositeMatrixAsTreeMatrix& other) { cm_ = other.m_->Clone(); }

  // Vector space of the Matrix's domain.
  virtual Teuchos::RCP<const TreeVectorSpace> domain() const
  {
    if (domain_ == Teuchos::null) domain_ = Teuchos::rcp(new TreeVectorSpace(cm_->domain()));
    return domain_;
  }

  // Vector space of the Matrix's range.
  virtual Teuchos::RCP<const TreeVectorSpace> range() const
  {
    if (range_ == Teuchos::null) range_ = Teuchos::rcp(new TreeVectorSpace(cm_->range()));
    return range_;
  }

  // Virtual copy constructor.
  virtual Teuchos::RCP<TreeMatrix> Clone() const
  {
    return Teuchos::rcp(new CompositeMatrixAsTreeMatrix(*this));
  }

  // Apply matrix, b <-- Ax, returns ierr
  virtual int Apply(const TreeVector& x, TreeVector& b) const
  {
    return cm_->Apply(x.data(), b.data());
  }

  // Apply the inverse, x <-- A^-1 b, returns ierr
  virtual int ApplyInverse(const TreeVector& b, TreeVector& x) const
  {
    return cm_->ApplyInverse(x.data(), b.data());
  }
};

} // namespace Amanzi
