/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for BlockVector, an implementation of a slightly improved
   Epetra_MultiVector which spans multiple simplices and knows how to
   communicate itself.
   ------------------------------------------------------------------------- */

#ifndef AMANZI_BLOCKVECTOR_HH_
#define AMANZI_BLOCKVECTOR_HH_

#include <vector>
#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"

#include "dbc.hh"
#include "data_structures_types.hh"

namespace Amanzi {

class BlockVector {

public:
  // Constructor
  BlockVector(const Epetra_MpiComm* comm,
                  std::vector<std::string>& names,
                  std::vector<Teuchos::RCP<const Epetra_Map> >& maps,
                  std::vector<int> num_dofs);

  // copy constructor
  BlockVector(const BlockVector& other, ConstructMode mode=CONSTRUCT_WITH_NEW_DATA);

  // assignment
  BlockVector& operator=(const BlockVector& other);

  // Check consistency of meta-data and allocate data.
  // Separating this from the object constructor is an active choice -- in
  // this way the non-CreateData version can represent just the structure of
  // the vector, much like an Epetra_Map.
  void CreateData();

  // Accessors

  // Accessors to meta-data
  // Iteration over names of the vector
  typedef std::vector<std::string>::const_iterator name_iterator;
  name_iterator begin() { return names_.begin(); }
  name_iterator end() { return names_.end(); }

  int num_components() const { return num_components_; }
  int num_dofs(std::string name) const { return num_dofs_[index_(name)]; }
  int size(std::string name) const { return sizes_[index_(name)]; }
  Teuchos::RCP<const Epetra_Map> map(std::string name) const { return maps_[index_(name)]; }


  // Accessors to data.
  // -- Access a view of a single component's data.
  Teuchos::RCP<const Epetra_MultiVector>
  ViewComponent(std::string name) const;

  Teuchos::RCP<Epetra_MultiVector>
  ViewComponent(std::string name);

  // -- View entries in the vectors.
  double operator()(std::string name, int i, int j) const {
    return (*data_[index_(name)])[i][j];
  }
  double operator()(std::string name, int j) const {
    return (*data_[index_(name)])[0][j];
  }

  // Mutators of data
  // -- Set entries in the vectors.
  void SetComponent(std::string name, const Teuchos::RCP<Epetra_MultiVector>& data);

  double& operator()(std::string name, int i, int j) {
    return (*data_[index_(name)])[i][j];
  }
  double& operator()(std::string name, int j) {
    return (*data_[index_(name)])[0][j];
  }

  // Vector operations.
  // Insert value into data.
  int PutScalar(double scalar);

  // Insert values into data, by DOF, not by component!
  int PutScalar(std::vector<double> scalar);

  // Insert value into component [name].
  int PutScalar(std::string name, double scalar);

  // Insert values into component [name].
  int PutScalar(std::string name, std::vector<double> scalar);

  // this <- this * scalarThis
  int Scale(double value);

  // Scale() applied to component name.
  int Scale(std::string name, double scalarThis);

  // this <- this + scalarA
  int Shift(double scalarA);

  // Shift() applied to component name.
  int Shift(std::string name, double scalarA);

  // result <- other \dot this
  int Dot(const BlockVector& other, double* result) const;

  // this <- scalarA*A + scalarThis*this
  BlockVector& Update(double scalarA, const BlockVector& A, double scalarThis);

  // this <- scalarA*A + scalarB*B + scalarThis*this
  BlockVector& Update(double scalarA, const BlockVector& A,
                          double scalarB, const BlockVector& B, double scalarThis);

  // this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
  int Multiply(double scalarAB, const BlockVector& A, const BlockVector& B,
               double scalarThis);

  // this <- scalarAB * B/A + scalarThis*this  (/ is the elementwise division
  int ReciprocalMultiply(double scalarAB, const BlockVector& A, const BlockVector& B,
                         double scalarThis);

  // Norms.
  int NormInf(double* norm) const;
  int Norm1(double* norm) const;
  int Norm2(double* norm) const;

  // Extras
  void Print(ostream &os) const;

private:
  int index_(std::string name) const {
    std::map<std::string, int>::const_iterator item = indexmap_.find(name);
    ASSERT(item != indexmap_.end());
    return item->second;
  }

private:
  const Epetra_MpiComm* comm_;
  std::map< std::string, int > indexmap_;
  int num_components_;

  std::vector<std::string> names_;
  std::vector<int> num_dofs_;
  std::vector<Teuchos::RCP<const Epetra_Map> > maps_;
  std::vector<int> sizes_;
  std::vector<Teuchos::RCP<Epetra_MultiVector> > data_;
};

} // namespace

#endif
