/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for TreeVector, a nested, hierarchical data structure for PK
   hiearchies.  This nested vector allows each physical PK to push back
   Epetra_MultiVectors to store their solution, and allows MPCs to push back
   TreeVectors in a nested format.  It also provides an implementation of the
   Vector interface for use with time integrators/nonlinear solvers.
   ------------------------------------------------------------------------- */

#ifndef DATA_STRUCTURES_TREEVECTOR_HH_
#define DATA_STRUCTURES_TREEVECTOR_HH_

#include <string>
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Vector.hh"

namespace Amanzi {

  class TreeVector : public Vector {

  public:
    // Basic constructor of an empty TreeVector
    explicit TreeVector(std::string name);

    // These constructors allocate their own memory and do NOT copy the values,
    // just the layout, of their arguments.

    // NOTE: are these really necessary?  I'm a little worried that the
    // semantics of these aren't clear and they likely won't be necessary.  I
    // think it's likely that TreeVectors get constructed in one of two ways:
    // via the standard copy constructor, or as an empty TreeVector that then
    // gets populated.  I suspect we should just remove these. --etc
    // TreeVector(std::string name, const Teuchos::RCP<Epetra_MultiVector>&);
    // TreeVector(std::string name, const Teuchos::RCP<Epetra_Vector>&);

    // these copy constructors are clearly necessary
    TreeVector(std::string name, const Teuchos::RCP<TreeVector>&);
    TreeVector(std::string name, const TreeVector&);
    TreeVector(const TreeVector&);
    TreeVector(const Teuchos::RCP<TreeVector>&);

    // set data
    //  - guaranteed to not error
    TreeVector& operator=(double value);

    //  - not guaranteed not to error
    TreeVector& operator=(const Epetra_Vector &value);
    TreeVector& operator=(const Epetra_MultiVector &value);
    TreeVector& operator=(const TreeVector &value);

    // metadata
    void SetName(std::string name) {
      name_ = name;
    }
    std::string Name() {
      return name_;
    }

    // operations
    // this <- scalar
    int PutScalar(double scalar);

    // ninf <- || this ||_{inf}
    int NormInf(double* ninf) const;

    // non-inherited extras
    void Print(ostream &os) const;

    // this <- value*this
    void Scale(double value);

    // this <- this + scalarA
    void Shift(double scalarA);

    // result <- other \dot this
    int Dot(const Vector& other, double* result) const;

    // this <- scalarA*A + scalarThis*this
    Vector& Update(double scalarA, const Vector& A, double scalarThis);

    // this <- scalarA*A + scalarB*B + scalarThis*this
    Vector& Update(double scalarA, const Vector& A,
                   double scalarB, const Vector& B, double scalarThis);

    // this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
    int Multiply(double scalarAB, const Vector& A, const Vector& B, double scalarThis);

    // Get a pointer to the sub-vector "subname".
    int SubVector(std::string subname, Teuchos::RCP<TreeVector>& subvec);
    int SubVector(std::string subname, Teuchos::RCP<const TreeVector>& subvec) const;

    // Get a pointer to the data vector indexed by vecnum.
    Teuchos::RCP<Epetra_MultiVector> operator[](int vecnum) {
      return data_[vecnum];
    }
    Teuchos::RCP<const Epetra_MultiVector> operator[](int vecnum) const {
      return data_[vecnum];
    }

    // Add a sub-vector as a child of this node.
    void PushBack(Teuchos::RCP<TreeVector>& subvec);

    // Add a data vector to this node of the tree.
    void PushBack(Teuchos::RCP<Epetra_MultiVector>& data);

  private:
    std::string name_;
    std::vector< Teuchos::RCP<Epetra_MultiVector> > data_;
    std::vector< Teuchos::RCP<TreeVector> > subvecs_;
  };

} // namespace


#endif
