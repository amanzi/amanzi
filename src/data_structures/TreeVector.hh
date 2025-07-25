/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/* -------------------------------------------------------------------------
ATS / Amanzi

Interface for TreeVector, a nested, hierarchical data structure for PK
hiearchies.  This vector allows each physical PK to use CompositeVectors to
store their solution, and allows MPCs to push back TreeVectors in a tree
format.

This vector provides the duck-type interface Vec and may be used with time
integrators/nonlinear solvers.

Note that a TreeVector may have EITHER subvecs OR data (in the form of a
CompositeVector), but NOT BOTH!  Nothing HERE precludes this, but it is
assumed in several places.
------------------------------------------------------------------------- */

#ifndef DATA_STRUCTURES_TREEVECTOR_HH_
#define DATA_STRUCTURES_TREEVECTOR_HH_

#include <string>
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

//#include "Iterators.hh"
#include "data_structures_types.hh"

#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"
#include "TreeVectorSpace.hh"

namespace Amanzi {

class TreeVector {
 public:
  using VectorSpace_t = TreeVectorSpace;
  // -- Constructors --

  // Basic constructors of a TreeVector
  TreeVector();
  explicit TreeVector(const Comm_ptr_type& comm);
  explicit TreeVector(const TreeVectorSpace& space, InitMode mode = init_mode_default);
  explicit TreeVector(const Teuchos::RCP<TreeVectorSpace>& space,
                      InitMode mode = init_mode_default);

  // copy constructors
  TreeVector(const TreeVector& other);

  // Assignment operator.
  TreeVector& operator=(const TreeVector& other);

  // -- Accessors --

  // Access to ANY communicator (this may be ill-posed!)
  Comm_ptr_type Comm() const { return Map().Comm(); }

  // Access to the space.
  const TreeVectorSpace& Map() const { return *map_; }
  const Teuchos::RCP<TreeVectorSpace>& get_map() const { return map_; }

  // Access to SubVectors
  using iterator = std::vector<Teuchos::RCP<TreeVector>>::iterator;
  using const_iterator = Teuchos::RCP<const TreeVector> const* const;

  // this is a very poor-man's hacky replacement for an iterator_adaptor
  // which aims to make const_iterators iterate over pointers to const objects.
  const_iterator begin() const
  {
    return reinterpret_cast<const Teuchos::RCP<const TreeVector>*>(&*subvecs_.begin());
  }
  const_iterator end() const
  {
    return reinterpret_cast<const Teuchos::RCP<const TreeVector>*>(&*subvecs_.end());
  }

  iterator begin() { return iterator(subvecs_.begin()); }
  iterator end() { return iterator(subvecs_.end()); }
  size_t size() const { return subvecs_.size(); }

  // Access to the sub-vector by index
  Teuchos::RCP<TreeVector> SubVector(int index);

  // Const access to the sub-vector by index
  Teuchos::RCP<const TreeVector> SubVector(int index) const;

  // Access to the data CompositeVector.
  Teuchos::RCP<CompositeVector> Data() { return data_; }

  // Const access to the data CompositeVector.
  Teuchos::RCP<const CompositeVector> Data() const { return data_; }


  // Add a sub-vector as a child of this node.
  void PushBack(const Teuchos::RCP<TreeVector>& subvec);

  // Set data by pointer
  void SetData(const Teuchos::RCP<CompositeVector>& data);


  // -- Assorted vector operations, this implements a Vec --
  // total length of the containing data
  int GlobalLength() const;

  // this <- scalar
  int PutScalar(double scalar);
  int PutScalarMasterAndGhosted(double scalar);
  int PutScalarGhosted(double scalar);

  // this <- random
  int Random();

  // n_l <- || this ||_{l}
  int Norm2(double* n2) const;
  int Norm1(double* n1) const;
  int NormInf(double* ninf) const;

  // this <- abs(this)
  int Abs(const TreeVector& other);

  // this <- value*this
  int Scale(double value);

  // this <- this + scalarA
  int Shift(double scalarA);

  // this <- element wise reciprocal(this)
  int Reciprocal(const TreeVector& other);

  // result <- other \dot this
  int Dot(const TreeVector& other, double* result) const;

  // this <- scalarA*A + scalarThis*this
  TreeVector& Update(double scalarA, const TreeVector& A, double scalarThis);

  // this <- scalarA*A + scalarB*B + scalarThis*this
  TreeVector& Update(double scalarA,
                     const TreeVector& A,
                     double scalarB,
                     const TreeVector& B,
                     double scalarThis);

  // this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
  int Multiply(double scalarAB, const TreeVector& A, const TreeVector& B, double scalarThis);

  // this <- scalarAB * A^-1@B + scalarThis*this  (@ is the elementwise product
  int ReciprocalMultiply(double scalarAB,
                         const TreeVector& A,
                         const TreeVector& B,
                         double scalarThis);

  // non-inherited extras
  void Print(std::ostream& os, bool data_io = true) const;

  // int GlobalLength() { std::cerr << "This method is not yet implemented\n"; return 0; }

 private:
  // Init's version of PushBack, which does not add to the space.
  void InitPushBack_(const Teuchos::RCP<TreeVector>& subvec);
  void InitMap_(InitMode mode);

 private:
  Teuchos::RCP<TreeVectorSpace> map_;

  Teuchos::RCP<CompositeVector> data_;
  std::vector<Teuchos::RCP<TreeVector>> subvecs_;
};


// non-member functions
inline Teuchos::RCP<TreeVector>
CreateTVwithOneLeaf(Teuchos::RCP<CompositeVector> cv)
{
  auto tvs = CreateTVSwithOneLeaf(cv->Map());
  auto tv = Teuchos::rcp(new TreeVector(tvs));
  tv->SetData(cv);
  return tv;
}

} // namespace Amanzi


#endif
