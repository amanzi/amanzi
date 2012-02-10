/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for CompositeVector, an implementation of a slightly improved
   Epetra_MultiVector which spans multiple simplices and knows how to
   communicate itself.

   NOTE: All CompositeVector data is NOT initialized to zero!
   ------------------------------------------------------------------------- */

#include <vector>
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CombineMode.h"

#include "Mesh.hh"
#include "MeshDefs.hh"
#include "CompositeVector.hh"

//namespace Amanzi {

CompositeVector::CompositeVector(Teuchos::RCP<AmanziMesh::Mesh>& mesh,
        std::vector<std::string> names,
        std::vector<AmanziMesh::Entity_kind> locations,
        std::vector<int> num_dofs, bool ghosted=true) :
    mesh_(mesh), names_(names), locations_(locations),
    num_dofs_(num_dofs), ghosted_(ghosted) {

  num_components_ = locations_.size();
  VerifyAndCreateData_();
};

CompositeVector::CompositeVector(Teuchos::RCP<AmanziMesh::Mesh>& mesh,
        std::vector<std::string> names,
        std::vector<AmanziMesh::Entity_kind> locations,
        int num_dofs, bool ghosted=true) :
    mesh_(mesh), names_(names), locations_(locations), ghosted_(ghosted) {

  num_components_ = locations_.size();
  num_dofs_.resize(num_components_);
  for (i = 0; i < num_components_; ++i) num_dofs_[i] = num_dofs;
  VerifyAndCreateData_();
};

CompositeVector::CompositeVector(Teuchos::RCP<AmanziMesh::Mesh>& mesh,
        std::string name, AmanziMesh::Entity_kind location,
        int num_dofs, bool ghosted=true) :
    mesh_(mesh), ghosted_(ghosted) {

  num_components_ = 1;
  num_dofs_.resize(1);
  num_dofs[0] = num_dofs;
  names_.resize(1);
  names_[0] = name;
  locations_.resize(1);
  locations_[0] = location;

  VerifyAndCreateData_();
};

CompositeVector::CompositeVector(const CompositeVector& other) :
    mesh_(other.mesh_), names_(other.names_),
    locations_(other.locations_), num_dofs_(other.num_dofs_),
    cardinalities_(other.cardinalities_), importers_(other.importers_),
    exporters_(other.exporters_), subfield_names_(other.subfield_names_),
    ghosted_(other.ghosted_), indexmap_(other.indexmap_) {

  num_components_ = other.num_components_;
  data_.resize(num_components_);
  owned_data_.resize(num_components_);

  for (unsigned int i = 0; i < num_components_; ++i) {
    data_[i] = Teuchos::rcp(new Epetra_MultiVector(*other.data_[i]));
  }

  // if the vector is not ghosted, we don't need owned_data_, as it is the
  // same as data_.
  if (!ghosted) {
    owned_data_ = data_;
  }
};

CompositeVector::CompositeVector(const CompositeVector& other, ConstructMode mode) :
    mesh_(other.mesh_), names_(other.names_), locations_(other.locations_),
    num_dofs_(other.num_dofs_), indexmap_(other.indexmap_),
    subfield_names_(other.subfield_names_, ghosted_(other.ghosted_) {

  num_components_ = other.num_components_;
  data_.resize(num_components_);
  owned_data_.resize(num_components_);
  cardinalities_.resize(num_components_);
  importers_.resize(num_components_);
  exporters_.resize(num_components_);
};

// Sets sizes of vectors, instantiates Epetra_Vectors, and preps for lazy
// creation of everything else.
void VerifyAndCreateData_() {
  ASSERT(names_.size() == num_components_);
  ASSERT(num_dofs_.size() == num_components_);

  // set sizes of vectors
  data_.resize(num_components_);
  owned_data_.resize(num_components_);
  cardinalities_.resize(num_components_);
  importers_.resize(num_components_);
  exporters_.resize(num_components_);
  subfield_names_.resize(num_components_);

  // create the map and vector data
  for (unsigned int i = 0; i < num_components_; ++i) {
    indexmap_[names_[i]] = i;

    if (locations_[i] == AmanziMesh::CELL) {
      Epetra_Map map = mesh_.cell_map(ghosted_);
      if (1 == num_dofs_[i]) {
        data_[i] = Teuchos::rcp(Epetra_Vector(map, false));
      } else {
        data_[i] = Teuchos::rcp(Epetra_MultiVector(map, num_dofs_[i], false));
      }
    } else if (locations_[i] == AmanziMesh::FACE) {
      Epetra_Map map = mesh_.face_map(ghosted_);
      if (1 == num_dofs_[i]) {
        data_[i] = Teuchos::rcp(Epetra_Vector(map, false));
      } else {
        data_[i] = Teuchos::rcp(Epetra_MultiVector(map, num_dofs_[i], false));
      }
    } else if (locations_[i] == AmanziMesh::NODE) {
      Epetra_Map map = mesh_.node_map(ghosted_);
      if (1 == num_dofs_[i]) {
        data_[i] = Teuchos::rcp(Epetra_Vector(map, false));
      } else {
        data_[i] = Teuchos::rcp(Epetra_MultiVector(map, num_dofs_[i], false));
      }
    } else {
      ASSERT(false);
    }

    cardinalities_[i] = map.NumMyElements();
    subfield_names_[i].resize(num_dofs_[i], std::string(""));
  }

  // if the vector is not ghosted, we don't need owned_data_, as it is the
  // same as data_.
  if (!ghosted_) {
    owned_data_ = data_;
  }
};


CompositeVector& CompositeVector::operator=(const CompositeVector& other) {
  if (this != &other) {
    if ((num_components_ != other.num_components_) || (mesh_ != other.mesh_) ||
        (ghosted_ != other.ghosted_)) {
      Errors::Message message("Attempted assignment of non-compatible CompositeVectors.");
      Exceptions::amanzi_throw(message);
    }
    for (unsigned int i = 0; i < num_components_; ++i) {
      if ((num_dofs_[i] != other.num_dofs_[i]) ||
          (cardinalities_[i] != other.cardinalities_[i])) {
        Errors::Message message("Attempted assignment of non-compatible CompositeVectors.");
        Exceptions::amanzi_throw(message);
      }
    }

    for (unsigned int i = 0; i < num_components_; ++i) {
      *data_[i] = *other.data_[i];
    }
    owned_composite_ = Teuchos::null;
  }
  return *this;
};

// view data
// -- Access a view of a single component's data.
// Ghosted views are simply the vector itself, while non-ghosted views are
// lazily generated.

// this one is broken without mutable?
Teuchos::RCP<const Epetra_MultiVector> CompositeVector::ViewComponent(std::string name,
        bool ghosted) const {
  unsigned int index;
  if (name == "" && num_components_ == 1) {
    index = 0;
  } else {
    index = indexmap_[name];
  }

  if (ghosted) {
    return data_[i];
  } else {
    if (owned_data_[i] == Teuchos::null) {
      CreateOwnedView_(i);
    }
    return owned_data_[i];
  }
};

Teuchos::RCP<Epetra_MultiVector> CompositeVector::ViewComponent(std::string name,
        bool ghosted) {
  unsigned int index;
  if (name == "" && num_components_ == 1) {
    index = 0;
  } else {
    index = indexmap_[name];
  }

  if (ghosted) {
    return data_[i];
  } else {
    if (owned_data_[i] == Teuchos::null) {
      CreateOwnedView_(i);
    }
    return owned_data_[i];
  }
};

// -- Access a view of the owned composite data.
// All meta-data is copied, but pointers to the data are replaced by
// pointers to owned data.

// this one is broken without mutable?
Teuchos::RCP<const CompositeVector> CompositeVector::ViewOwned() const {
  if (owned_composite_ == Teuchos::null) {
    // create the owned view
    owned_composite_ = Teuchos::rcp(new CompositeVector(*this, CONSTRUCT_WITHOUT_DATA));
    owned_composite_->ghosted_ = false;

    if (ghosted_) {
      for (unsigned int i = 0; i < num_components_; ++i) {
        owned_composite_->data_[i] = ViewComponent(names_[i], false);
      }
    } else {
      owned_composite_->data_ = data_;
    }
    owned_composite_->owned_data_ = owned_composite_->data_;
  }
  return owned_composite_;
};

Teuchos::RCP<CompositeVector> CompositeVector::ViewOwned() {
  if (owned_composite_ == Teuchos::null) {
    // create the owned view
    owned_composite_ = Teuchos::rcp(new CompositeVector(*this, CONSTRUCT_WITHOUT_DATA));
    owned_composite_->ghosted_ = false;

    if (ghosted_) {
      for (unsigned int i = 0; i < num_components_; ++i) {
        owned_composite_->data_[i] = ViewComponent(names_[i], false);
      }
    } else {
      owned_composite_->data_ = data_;
    }
    owned_composite_->owned_data_ = owned_composite_->data_;
  }
  return owned_composite_;
};


  ///////// DONE THROUGH HERE /////////


  // communicate
    // -- Scatter master values to ghosted values.
    // Modes shown in Epetra_CombineMode.h, but the default is Insert, which
    // overwrites the current ghost value with the (unique) new master value.
    void CompositeVector::ScatterMasterToGhosted(Epetra_CombineMode mode=Insert);

    // -- Combine ghosted values back to master values.
    // Modes shown in Epetra_CombineMode.h, but the default is InsertAdd,
    // where off-process values are first summed, then replace the current
    // value.
    void CompositeVector::GatherGhostedToMaster(Epetra_CombineMode mode=InsertAdd);

    // assorted vector operations for use by time integrators?
    // These may make life much easier for time integration of the full
    // CompositeVector.  For instance, one could implement the BDF integrators
    // doing all operations by first calling ViewOwned() to get a non-ghosted
    // CompositeVector, and then using these operations.

    // -- Insert value into data.
    int CompositeVector::PutScalar(double scalar);
    int CompositeVector::PutScalar(std::vector<double> scalar);
    int CompositeVector::PutScalar(int index, double scalar);
    int CompositeVector::PutScalar(int index, std::vector<double> scalar);
    int CompositeVector::PutScalar(int index, int blockid, double scalar);
    int CompositeVector::PutScalar(int index, int blockid, std::vector<double> scalar);

    // -- this <- value*this
    void CompositeVector::Scale(double value);

    // -- this <- this + scalarA
    void CompositeVector::Shift(double scalarA);

    // -- result <- other \dot this
    int CompositeVector::Dot(const CompositeVector& other, double* result) const;

    // -- this <- scalarA*A + scalarThis*this
    CompositeVector& CompositeVector::Update(double scalarA, const CompositeVector& A, double scalarThis);

    // -- this <- scalarA*A + scalarB*B + scalarThis*this
    CompositeVector& CompositeVector::Update(double scalarA, const CompositeVector& A,
                            double scalarB, const CompositeVector& B, double scalarThis);

    // -- this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
    int CompositeVector::Multiply(double scalarAB, const CompositeVector& A, const CompositeVector& B,
                 double scalarThis);

    // -- norms
    int CompositeVector::NormInf(double* norm) const;
    int CompositeVector::NormTwo(double* norm) const;

  };

} // namespace

#endif
