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

*/

#ifndef DATA_STRUCTURES_TREEVECTOR_HH_
#define DATA_STRUCTURES_TREEVECTOR_HH_

#include <string>
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Teuchos_DataAccess.hpp"

#include "DataStructuresHelpers.hh"

#include "CompositeVector.hh"
#include "TreeVectorSpace.hh"

namespace Amanzi {

class TreeVector {
 public:
  using VectorSpace_t = TreeVectorSpace;
  // -- Constructors --
  // Basic constructors of a TreeVector
  explicit TreeVector(const Teuchos::RCP<const TreeVectorSpace>& space,
                      InitMode mode = InitMode::ZERO);

  // copy constructors
  TreeVector(const TreeVector& other,
             Teuchos::DataAccess access = Teuchos::DataAccess::Copy,
             InitMode mode = InitMode::COPY);

  // Assignment operator.
  TreeVector& operator=(const TreeVector& other);
  void assign(const TreeVector& other) { *this = other; }

  // -- Accessors --

  // Access to ANY communicator (this may be ill-posed!)
  Comm_ptr_type getComm() const { return getMap()->getComm(); }

  // Access to the space.
  const Teuchos::RCP<const TreeVectorSpace>& getMap() const { return map_; }

  // Access to SubVectors
  using iterator = std::vector<Teuchos::RCP<TreeVector>>::iterator;
  using const_iterator = Teuchos::RCP<const TreeVector> const * const;

  // this is a very poor-man's hacky replacement for boost iterator_adaptor
  // which aims to make const_iterators iterate over pointers to const objects.
  const_iterator begin() const {
    return reinterpret_cast<const Teuchos::RCP<const TreeVector>*>(&*subvecs_.begin());
  }
  const_iterator end() const {
    return reinterpret_cast<const Teuchos::RCP<const TreeVector>*>(&*subvecs_.end());
  }

  iterator begin() { return iterator(subvecs_.begin()); }
  iterator end() { return iterator(subvecs_.end()); }
  size_t size() const { return subvecs_.size(); }

  // Access to the sub-vector by index
  Teuchos::RCP<TreeVector> getSubVector(int index);

  // Const access to the sub-vector by index
  Teuchos::RCP<const TreeVector> getSubVector(int index) const;

  // Access to the sub-vector by index
  void setSubVector(int index, const Teuchos::RCP<TreeVector>& tv);

  // Access to the data CompositeVector.
  Teuchos::RCP<CompositeVector> getData() { return data_; }

  // Const access to the data CompositeVector.
  Teuchos::RCP<const CompositeVector> getData() const { return data_; }

  // Set data by pointer
  void setData(const Teuchos::RCP<CompositeVector>& data);


  // -- Assorted vector operations, this implements a Vec --
  // total length of the containing data
  int getGlobalLength() const;

  // this <- scalar
  void putScalar(double scalar);
  void putScalarMasterAndGhosted(double scalar);
  void putScalarGhosted(double scalar);

  // this <- random
  void randomize();

  // n_l <- || this ||_{l}
  double norm2() const;
  double norm1() const;
  double normInf() const;

  // this <- abs(this)
  void abs(const TreeVector& other);

  // this <- value*this
  void scale(double value);

  // this <- this + scalarA
  void shift(double scalarA);

  // this <- element wise reciprocal(this)
  void reciprocal(const TreeVector& other);

  // result <- other \dot this
  double dot(const TreeVector& other) const;

  // this <- scalarA*A + scalarThis*this
  void update(double scalarA, const TreeVector& A, double scalarThis);

  // this <- scalarA*A + scalarB*B + scalarThis*this
  void update(double scalarA,
              const TreeVector& A,
              double scalarB,
              const TreeVector& B,
              double scalarThis);

  // this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
  void
  elementWiseMultiply(double scalarAB, const TreeVector& A, const TreeVector& B, double scalarThis);

  // this <- scalarAB * A^-1@B + scalarThis*this  (@ is the elementwise product
  // int ReciprocalelementWiseMultiply(double scalarAB, const TreeVector& A,
  //                        const TreeVector& B, double scalarThis);

  // non-inherited extras
  void print(std::ostream& os) const;

  int getGlobalLength()
  {
    std::cerr << "This method is not yet implemented\n";
    return 0;
  }

 private:
  // Init's version of PushBack, which does not add to the space.
  void InitPushBack_(const Teuchos::RCP<TreeVector>& subvec);
  void InitMap_(InitMode mode);

 private:
  Teuchos::RCP<const TreeVectorSpace> map_;

  Teuchos::RCP<CompositeVector> data_;
  std::vector<Teuchos::RCP<TreeVector>> subvecs_;
};

} // namespace Amanzi


#endif
