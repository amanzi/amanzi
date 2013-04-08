/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for CompositeVector, an implementation of a slightly improved
   Epetra_MultiVector which spans multiple simplices and knows how to
   communicate itself.
   ------------------------------------------------------------------------- */

#ifndef COMPOSITEVECTOR_HH_
#define COMPOSITEVECTOR_HH_

#include <vector>
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_CombineMode.h"
#include "Epetra_Import.h"

#include "dbc.hh"
#include "Mesh.hh"
#include "data_structures_types.hh"
#include "block_vector.hh"

namespace Amanzi {

class CompositeVector {

public:
  // Constructors
  CompositeVector(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      mesh_(mesh), created_(false), ghosted_(false) {}

  CompositeVector(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                  std::vector<std::string> names,
                  std::vector<AmanziMesh::Entity_kind> locations,
                  std::vector<int> num_dofs,
                  bool ghosted);

  // copy constructor / assignment
  CompositeVector(const CompositeVector& other,
                  ConstructMode mode=CONSTRUCT_WITH_NEW_DATA);

  CompositeVector& operator=(const CompositeVector& other);

  // Check consistency of meta-data and allocate data.
  virtual void CreateData();

  // Accessors to meta-data
  // -- Iteration over names of the vector
  typedef std::vector<std::string>::const_iterator name_iterator;
  name_iterator begin() const { return names_.begin(); }
  name_iterator end() const { return names_.end(); }

  bool has_component(std::string name) const {
    return indexmap_.find(name) != indexmap_.end(); }

  int num_components() const { return num_components_; }
  Teuchos::RCP<const AmanziMesh::Mesh> mesh() const { return mesh_; }
  const Epetra_MpiComm* comm() { return mesh_->get_comm(); }

  int num_dofs(std::string name) const { return ghostvec_->num_dofs(name); }
  AmanziMesh::Entity_kind location(std::string name) const {
    return locations_[index_(name)];
  }

  bool ghosted() const { return ghosted_; }
  bool created() const { return created_; }

  int size(std::string name, bool ghosted=false) const {
    return ghosted ? ghostvec_->size(name) : mastervec_->size(name); }

  // View data.
  // -- Access the full block system -- likely not often needed.
  Teuchos::RCP<const BlockVector> get_master_vector() const { return mastervec_; }
  Teuchos::RCP<BlockVector> get_master_vector() { return mastervec_; }
  Teuchos::RCP<const BlockVector> get_ghost_vector() const { return ghostvec_; }
  Teuchos::RCP<BlockVector> get_ghost_vector() { return ghostvec_; }

  // -- Access a view of a single component's data.
  Teuchos::RCP<const Epetra_Map> map(std::string name, bool ghosted=false) const {
    return ghosted ? ghostvec_->map(name) : mastervec_->map(name);
  }

  Teuchos::RCP<const Epetra_MultiVector>
  ViewComponent(std::string name, bool ghosted=false) const;

  Teuchos::RCP<Epetra_MultiVector>
  ViewComponent(std::string name, bool ghosted=false);

  // -- view entries in the vectors
  double operator()(std::string name, int i, int j) const {
    return (*ghostvec_)(name,i,j);
  }
  double operator()(std::string name, int j) const {
    return (*ghostvec_)(name,0,j);
  }

  // Set data.
  // -- Set block by pointer if possible, copy if not.
  void SetComponent(std::string name, const Teuchos::RCP<Epetra_MultiVector>& data);

  // -- set entries in the vectors
  double& operator()(std::string name, int i, int j) {
    return (*ghostvec_)(name,i,j);
  }
  double& operator()(std::string name, int j) {
    return (*ghostvec_)(name,0,j);
  }

  // Communicate data.
  // -- Scatter master values to ghosted values.
  // Modes shown in Epetra_CombineMode.h, but the default is Insert, which
  // overwrites the current ghost value with the (unique) new master value.
  //
  // Note that although scatter changes things, it doesn't change master
  // data, so we allow it to work on const.  This is necessary for a
  // non-owning PK to communicate a non-owned vector.
  virtual void ScatterMasterToGhosted(Epetra_CombineMode mode=Insert) const;
  virtual void ScatterMasterToGhosted(std::string name, Epetra_CombineMode mode=Insert) const;

  // -- Combine ghosted values back to master values.
  // Modes shown in Epetra_CombineMode.h, but the default is Add,
  // where off-process values are first summed into the on-process value.
  virtual void GatherGhostedToMaster(Epetra_CombineMode mode=Add);
  virtual void GatherGhostedToMaster(std::string name, Epetra_CombineMode mode=Add);

  // Assorted vector operations.
  //   These may make life much easier for time integration of the full
  //   CompositeVector.  For instance, one could implement the BDF integrators
  //   doing all operations by first calling ViewOwned() to get a non-ghosted
  //   CompositeVector, and then using these operations.

  // -- Insert scalar into data.
  int PutScalar(double scalar) { return mastervec_->PutScalar(scalar); }
  int PutScalar(std::vector<double> scalar) { return mastervec_->PutScalar(scalar); }
  int PutScalar(std::string name, double scalar) {
    return mastervec_->PutScalar(name, scalar); }
  int PutScalar(std::string name, std::vector<double> scalar) {
    return mastervec_->PutScalar(name, scalar); }

  // -- this <- scalar*this
  int Scale(double scalar) { return mastervec_->Scale(scalar); }
  int Scale(std::string name, double scalar) { return mastervec_->Scale(name, scalar); }

  // -- this <- this + scalarA
  int Shift(double scalar) { return mastervec_->Shift(scalar); }
  int Shift(std::string name, double scalar) { return mastervec_->Shift(name, scalar); }

  // -- result <- other \dot this
  int Dot(const CompositeVector& other, double* result) const {
    return mastervec_->Dot(*other.get_master_vector(), result); }

  // -- this <- scalarA*A + scalarThis*this
  CompositeVector& Update(double scalarA, const CompositeVector& A, double scalarThis) {
    mastervec_->Update(scalarA, *A.get_master_vector(), scalarThis);
    return *this;
  }

  // -- this <- scalarA*A + scalarB*B + scalarThis*this
  CompositeVector& Update(double scalarA, const CompositeVector& A,
                          double scalarB, const CompositeVector& B, double scalarThis) {
    mastervec_->Update(scalarA, *A.get_master_vector(),
                       scalarB, *B.get_master_vector(), scalarThis);
    return *this;
  }


  // -- this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
  int Multiply(double scalarAB, const CompositeVector& A, const CompositeVector& B,
               double scalarThis) {
    return mastervec_->Multiply(scalarAB, *A.get_master_vector(), *B.get_master_vector(),
            scalarThis);
  }

  // -- this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
  int ReciprocalMultiply(double scalarAB, const CompositeVector& A, const CompositeVector& B,
               double scalarThis) {
    return mastervec_->ReciprocalMultiply(scalarAB, *A.get_master_vector(), *B.get_master_vector(),
            scalarThis);
  }


  // -- norms
  int NormInf(double* norm) const { return mastervec_->NormInf(norm); }
  int Norm1(double* norm) const { return mastervec_->Norm1(norm); }
  int Norm2(double* norm) const { return mastervec_->Norm2(norm); }

  // Extras
  void Print(ostream &os) const { return mastervec_->Print(os); }

 protected:
  int index_(std::string name) const {
    std::map<std::string, int>::const_iterator item = indexmap_.find(name);
    ASSERT(item != indexmap_.end());
    return item->second;
  }

  void AssertCreatedOrDie_() const;
  //  void AssertCommonStructureOrDie_(const CompositeVector& other) const;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  bool ghosted_;
  bool created_; // has data been allocated already

  // data enumerating the blocks
  int num_components_;
  std::map< std::string, int > indexmap_;
  std::vector<std::string> names_;
  std::vector<AmanziMesh::Entity_kind> locations_;

  // data containing the blocks
  mutable Teuchos::RCP<BlockVector> ghostvec_;
  Teuchos::RCP<BlockVector> mastervec_;

  mutable std::vector<Teuchos::RCP<Epetra_Import> > importers_;
};

} // namespace

#endif
