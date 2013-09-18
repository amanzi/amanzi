/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon (ecoon@lanl.gov)

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

#include "data_structures_types.hh"
#include "composite_vector.hh"

namespace Amanzi {

class TreeVector {

 public:
  // -- Constructors --

  // Basic constructor of an empty TreeVector
  TreeVector();

  // Basic constructor of an empty TreeVector with a name.
  explicit TreeVector(std::string name);

  // Copy constructor.
  //
  // Available modes are:
  //  CONSTRUCT_WITH_NEW_DATA : creates the vector with new uninitialized data
  //  CONSTRUCT_WITH_OLD_DATA : creates a new vector shell with pointers to
  //                            the same old data
  TreeVector(const TreeVector& other,
             ConstructMode mode=CONSTRUCT_WITH_NEW_DATA);

  // Copy constructor with a new name.
  //
  // Available modes are:
  //  CONSTRUCT_WITH_NEW_DATA : creates the vector with new uninitialized data
  //  CONSTRUCT_WITH_OLD_DATA : creates a new vector shell with pointers to
  //                            the same old data
  TreeVector(std::string name, const TreeVector& other,
             ConstructMode mode=CONSTRUCT_WITH_NEW_DATA);

  // Assignment operator.
  TreeVector& operator=(const TreeVector& other);

  // -- Accessors --

  // name accessor
  std::string name() { return name_; }

  // Access to ANY communicator (this may be ill-posed!)
  const Epetra_MpiComm& Comm() {
    if (data_ != Teuchos::null) return data_->Comm();
    if (subvecs_.size() > 0) return subvecs_[0]->Comm();
    ASSERT(0);
  }

  // Access to SubVectors
  std::vector< Teuchos::RCP<TreeVector> > SubVectors() { return subvecs_; }

  // Const access to SubVectors.
  const std::vector< Teuchos::RCP<TreeVector> > SubVectors() const {
    return subvecs_; }

  // Access to a sub-vector by name.
  Teuchos::RCP<TreeVector> SubVector(std::string subname);

  // Const access to a sub-vector by name.
  Teuchos::RCP<const TreeVector> SubVector(std::string subname) const;

  // Access to the sub-vector by index
  Teuchos::RCP<TreeVector> SubVector(int index);

  // Const access to the sub-vector by index
  Teuchos::RCP<const TreeVector> SubVector(int index) const;

  // Access to the data CompositeVector.
  Teuchos::RCP<CompositeVector> data() { return data_; }

  // Const access to the data CompositeVector.
  Teuchos::RCP<const CompositeVector> data() const { return data_; }

  // -- Mutators --

  // name mutator
  void set_name(std::string name) { name_ = name; }

  // Add a sub-vector as a child of this node.
  void PushBack(const Teuchos::RCP<TreeVector>& subvec);

  // Set data by pointer
  void set_data(const Teuchos::RCP<CompositeVector>& data) { data_ = data; }

  // Set data by copy
  void set_data(const CompositeVector& data) { *data_ = data; }

  // -- Assorted vector operations, this implements a Vec --

  // this <- scalar
  int PutScalar(double scalar);

  // n_l <- || this ||_{l}
  int Norm2(double* n2) const;
  int Norm1(double* n1) const;
  int NormInf(double* ninf) const;

  // this <- value*this
  int Scale(double value);

  // this <- this + scalarA
  int Shift(double scalarA);

  // result <- other \dot this
  int Dot(const TreeVector& other, double* result) const;

  // this <- scalarA*A + scalarThis*this
  TreeVector& Update(double scalarA, const TreeVector& A, double scalarThis);

  // this <- scalarA*A + scalarB*B + scalarThis*this
  TreeVector& Update(double scalarA, const TreeVector& A,
                     double scalarB, const TreeVector& B, double scalarThis);

  // this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
  int Multiply(double scalarAB, const TreeVector& A, const TreeVector& B,
               double scalarThis);

  // non-inherited extras
  void Print(ostream &os) const;


  private:
    std::string name_;
    Teuchos::RCP<CompositeVector> data_;
    std::vector< Teuchos::RCP<TreeVector> > subvecs_;
  };

} // namespace


#endif
