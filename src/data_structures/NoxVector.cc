/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: amklinv
*/

/*
 * NoxVector.cc
 *
 *  Created on: May 3, 2016
 */

#include <NoxVector.hh>
#include <TreeVector.hh>

namespace Amanzi {
template <class VectorClass>
NoxVector<VectorClass>::NoxVector(const NoxVector& other, NOX::CopyType type)
{
  switch (type) {
  case NOX::DeepCopy:
    vec_ = Teuchos::rcp(new VectorClass(*other.vec_, INIT_MODE_COPY));
    break;
  case NOX::ShapeCopy:
    vec_ = Teuchos::rcp(new VectorClass(*other.vec_, INIT_MODE_NONE));
    break;
  default:
    vec_ = Teuchos::rcp(new VectorClass(*other.vec_, INIT_MODE_NONE));
    break;
  }
}

template <>
NoxVector<CompositeVector>::NoxVector(const NoxVector& other, NOX::CopyType type)
{
  switch (type) {
  case NOX::DeepCopy:
    vec_ = Teuchos::rcp(new CompositeVector(*other.vec_, INIT_MODE_COPY));
    break;
  case NOX::ShapeCopy:
    vec_ = Teuchos::rcp(new CompositeVector(*other.vec_, INIT_MODE_NONE));
    break;
  default:
    vec_ = Teuchos::rcp(new CompositeVector(*other.vec_, INIT_MODE_NONE));
    break;
  }
}

template <>
NoxVector<TreeVector>::NoxVector(const NoxVector& other, NOX::CopyType type)
{
  switch (type) {
  case NOX::DeepCopy:
    vec_ = Teuchos::rcp(new TreeVector(*other.vec_, INIT_MODE_COPY));
    break;
  case NOX::ShapeCopy:
    vec_ = Teuchos::rcp(new TreeVector(*other.vec_, INIT_MODE_NONE));
    break;
  default:
    vec_ = Teuchos::rcp(new TreeVector(*other.vec_, INIT_MODE_NONE));
    break;
  }
}

template <>
NoxVector<Epetra_Vector>::NoxVector(const NoxVector& other, NOX::CopyType type)
{
  vec_ = Teuchos::rcp(new Epetra_Vector(*other.vec_));
}
} // namespace Amanzi
