/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Data Structures

  Interface for BlockVector, an implementation of a slightly improved
  Epetra_MultiVector which spans multiple simplices and knows how to
  communicate itself.
*/

#ifndef AMANZI_BLOCKVECTOR_HH_
#define AMANZI_BLOCKVECTOR_HH_

#include <vector>
#include "Teuchos_RCP.hpp"
#include "Epetra_BlockMap.h"
#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"

#include "AmanziComm.hh"
#include "dbc.hh"
#include "data_structures_types.hh"

namespace Amanzi {

class BlockVector {
 public:
  // Constructor
  BlockVector(const Comm_ptr_type& comm,
              const std::vector<std::string>& names,
              const std::vector<Teuchos::RCP<const Epetra_BlockMap>>& maps,
              const std::vector<int>& num_dofs);

  // copy constructor
  BlockVector(const BlockVector& other);

  // Constructor just does maps, this creates data.
  void CreateData();

  // assignment
  BlockVector& operator=(const BlockVector& other);

  // Accessors

  // Accessors to meta-data
  // Iteration over names of the vector
  typedef std::vector<std::string>::const_iterator name_iterator;
  name_iterator begin() const { return names_.begin(); }
  name_iterator end() const { return names_.end(); }
  unsigned int size() const { return names_.size(); }

  int GlobalLength() const;
  long int GetLocalElementCount() const;
  int NumComponents() const { return num_components_; }
  int NumVectors(const std::string& name) const { return num_dofs_[Index_(name)]; }
  unsigned int size(const std::string& name) const { return sizes_[Index_(name)]; }

  Teuchos::RCP<const Epetra_BlockMap> ComponentMap(const std::string& name) const
  {
    return maps_[Index_(name)];
  }

  // Accessors to data.
  bool HasComponent(const std::string& name) const { return indexmap_.find(name) != indexmap_.end(); }

  // -- Access a view of a single component's data.
  Teuchos::RCP<const Epetra_MultiVector> ViewComponent(const std::string& name) const;

  Teuchos::RCP<Epetra_MultiVector> ViewComponent(const std::string& name);

  // -- View entries in the vectors.
  double operator()(const std::string& name, int i, int j) const { return (*data_[Index_(name)])[i][j]; }
  double operator()(const std::string& name, int j) const { return (*data_[Index_(name)])[0][j]; }

  // Mutators of data
  // -- Set entries in the vectors.
  void SetComponent(const std::string& name, const Teuchos::RCP<Epetra_MultiVector>& data);

  // double& operator()(std::string name, int i, int j) {
  //   return (*data_[Index_(name)])[i][j];
  // }
  // double& operator()(std::string name, int j) {
  //   return (*data_[Index_(name)])[0][j];
  // }

  // Vector operations.
  // Insert value into data.
  int PutScalar(double scalar);

  // Insert values into data, by DOF, not by component!
  int PutScalar(const std::vector<double>& scalar);

  // Insert value into component [name].
  int PutScalar(const std::string& name, double scalar);

  // Insert values into component [name].
  int PutScalar(const std::string& name, const std::vector<double>& scalar);

  // this <- abs(this)
  int Abs(const BlockVector& other);

  // this <- this * scalarThis
  int Scale(double value);

  // Scale() applied to component name.
  int Scale(const std::string& name, double scalarThis);

  // this <- this + scalarA
  int Shift(double scalarA);

  // Shift() applied to component name.
  int Shift(const std::string& name, double scalarA);

  // this <- element wise reciprocal(this)
  int Reciprocal(const BlockVector& other);

  // result <- other \dot this
  int Dot(const BlockVector& other, double* result) const;

  // this <- scalarA*A + scalarThis*this
  BlockVector& Update(double scalarA, const BlockVector& A, double scalarThis);

  // this <- scalarA*A + scalarB*B + scalarThis*this
  BlockVector& Update(double scalarA,
                      const BlockVector& A,
                      double scalarB,
                      const BlockVector& B,
                      double scalarThis);

  // this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
  int Multiply(double scalarAB, const BlockVector& A, const BlockVector& B, double scalarThis);

  // this <- scalarAB * B/A + scalarThis*this  (/ is the elementwise division
  int ReciprocalMultiply(double scalarAB,
                         const BlockVector& A,
                         const BlockVector& B,
                         double scalarThis);

  // Norms.
  int NormInf(double* norm) const;
  int Norm1(double* norm) const;
  int Norm2(double* norm) const;

  int MinValue(double* value) const;
  int MaxValue(double* value) const;
  int MeanValue(double* value) const;

  // Extras
  void Print(std::ostream& os, bool data_io = true) const;

  int Random();

  Comm_ptr_type Comm() const { return comm_; }

 private:
  int Index_(const std::string& name) const
  {
    std::map<std::string, int>::const_iterator item = indexmap_.find(name);
    AMANZI_ASSERT(item != indexmap_.end());
    return item->second;
  }

 private:
  Comm_ptr_type comm_;
  std::map<std::string, int> indexmap_;
  int num_components_;

  std::vector<std::string> names_;
  std::vector<int> num_dofs_;
  std::vector<Teuchos::RCP<const Epetra_BlockMap>> maps_;
  std::vector<unsigned int> sizes_;
  std::vector<Teuchos::RCP<Epetra_MultiVector>> data_;
};

} // namespace Amanzi

#endif
