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

#include "dbc.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "MeshDefs.hh"
#include "CompositeVector.hh"

namespace Amanzi {

CompositeVector::CompositeVector(Teuchos::RCP<AmanziMesh::Mesh>& mesh,
        std::vector<std::string> names,
        std::vector<AmanziMesh::Entity_kind> locations,
        std::vector<int> num_dofs, bool ghosted) :
    mesh_(mesh), names_(names), locations_(locations),
    num_dofs_(num_dofs), ghosted_(ghosted) {

  num_components_ = locations_.size();
  VerifyAndCreateData_();
};

CompositeVector::CompositeVector(Teuchos::RCP<AmanziMesh::Mesh>& mesh,
        std::vector<std::string> names,
        std::vector<AmanziMesh::Entity_kind> locations,
        int num_dofs, bool ghosted) :
    mesh_(mesh), names_(names), locations_(locations), ghosted_(ghosted) {

  num_components_ = locations_.size();
  num_dofs_.resize(num_components_);
  for (unsigned int i = 0; i != num_components_; ++i) num_dofs_[i] = num_dofs;
  VerifyAndCreateData_();

  if (num_dofs == 1) {
    subfield_names_.resize(num_components_);
    for (unsigned int i = 0; i != num_components_; ++i) {
      subfield_names_[i].resize(1);
      subfield_names_[i][0] = names[i];
    }
  }
};

CompositeVector::CompositeVector(Teuchos::RCP<AmanziMesh::Mesh>& mesh,
        std::string name, AmanziMesh::Entity_kind location,
        int num_dofs, bool ghosted) :
    mesh_(mesh), ghosted_(ghosted) {

  num_components_ = 1;
  num_dofs_.resize(1);
  num_dofs_[0] = num_dofs;
  names_.resize(1);
  names_[0] = name;
  locations_.resize(1);
  locations_[0] = location;

  VerifyAndCreateData_();

  if (num_dofs == 1) {
    subfield_names_.resize(num_components_);
    subfield_names_[0].resize(1);
    subfield_names_[0][0] = name;
  }
};

CompositeVector::CompositeVector(const CompositeVector& other) :
    mesh_(other.mesh_), names_(other.names_), locations_(other.locations_),
    num_dofs_(other.num_dofs_), cardinalities_(other.cardinalities_), importers_(other.importers_),
    subfield_names_(other.subfield_names_), ghosted_(other.ghosted_), indexmap_(other.indexmap_) {

  num_components_ = other.num_components_;
  data_.resize(num_components_);
  owned_data_.resize(num_components_);

  for (unsigned int i = 0; i != num_components_; ++i) {
    data_[i] = Teuchos::rcp(new Epetra_MultiVector(*other.data_[i]));
  }

  // if the vector is not ghosted, we don't need owned_data_, as it is the
  // same as data_.
  if (!ghosted_) {
    owned_data_ = data_;
  }
};

CompositeVector::CompositeVector(const CompositeVector& other, ConstructMode mode) :
    mesh_(other.mesh_), names_(other.names_), locations_(other.locations_),
    num_dofs_(other.num_dofs_), indexmap_(other.indexmap_),
    subfield_names_(other.subfield_names_), ghosted_(other.ghosted_) {

  num_components_ = other.num_components_;
  data_.resize(num_components_);
  owned_data_.resize(num_components_);
  cardinalities_.resize(num_components_);
  importers_.resize(num_components_);
};

// Sets sizes of vectors, instantiates Epetra_Vectors, and preps for lazy
// creation of everything else.
void CompositeVector::VerifyAndCreateData_() {
  ASSERT(names_.size() == num_components_);
  ASSERT(num_dofs_.size() == num_components_);

  // set sizes of vectors
  data_.resize(num_components_);
  owned_data_.resize(num_components_);
  cardinalities_.resize(num_components_);
  importers_.resize(num_components_);
  subfield_names_.resize(num_components_);

  // create the map and vector data
  for (unsigned int i = 0; i != num_components_; ++i) {
    indexmap_[names_[i]] = i;
    if (locations_[i] == AmanziMesh::CELL) {
      Epetra_Map map = mesh_->cell_map(ghosted_);
      cardinalities_[i] = map.NumMyElements();
      if (1 == num_dofs_[i]) {
        data_[i] = Teuchos::rcp(new Epetra_Vector(map, false));
      } else {
        data_[i] = Teuchos::rcp(new Epetra_MultiVector(map, num_dofs_[i], false));
      }
    } else if (locations_[i] == AmanziMesh::FACE) {
      Epetra_Map map = mesh_->face_map(ghosted_);
      cardinalities_[i] = map.NumMyElements();
      if (1 == num_dofs_[i]) {
        data_[i] = Teuchos::rcp(new Epetra_Vector(map, false));
      } else {
        data_[i] = Teuchos::rcp(new Epetra_MultiVector(map, num_dofs_[i], false));
      }
    } else if (locations_[i] == AmanziMesh::NODE) {
      Epetra_Map map = mesh_->node_map(ghosted_);
      cardinalities_[i] = map.NumMyElements();
      if (1 == num_dofs_[i]) {
        data_[i] = Teuchos::rcp(new Epetra_Vector(map, false));
      } else {
        data_[i] = Teuchos::rcp(new Epetra_MultiVector(map, num_dofs_[i], false));
      }
    } else {
      ASSERT(false);
    }

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
    for (unsigned int i = 0; i != num_components_; ++i) {
      if ((num_dofs_[i] != other.num_dofs_[i]) ||
          (cardinalities_[i] != other.cardinalities_[i])) {
        Errors::Message message("Attempted assignment of non-compatible CompositeVectors.");
        Exceptions::amanzi_throw(message);
      }
    }

    for (unsigned int i = 0; i != num_components_; ++i) {
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

// create the non-ghosted view
void CompositeVector::CreateOwnedView_(unsigned int i) const {
  double** data;
  data_[i]->ExtractView(&data);

  if (locations_[i] == AmanziMesh::CELL) {
    const Epetra_BlockMap& source_map = mesh_->cell_map(false);
    owned_data_[i] = Teuchos::rcp(new Epetra_MultiVector(View, source_map, data, num_dofs_[i]));
  } else if (locations_[i] == AmanziMesh::FACE) {
    const Epetra_BlockMap& source_map = mesh_->face_map(false);
    owned_data_[i] = Teuchos::rcp(new Epetra_MultiVector(View, source_map, data, num_dofs_[i]));
  } else if (locations_[i] == AmanziMesh::NODE) {
    const Epetra_BlockMap& source_map = mesh_->node_map(false);
    owned_data_[i] = Teuchos::rcp(new Epetra_MultiVector(View, source_map, data, num_dofs_[i]));
  } else {
    ASSERT(false);
  }
};

// this one is broken without mutable?
Teuchos::RCP<const Epetra_MultiVector> CompositeVector::ViewComponent(std::string name,
        bool ghosted) const {
  unsigned int i;
  if (name == "" && num_components_ == 1) {
    i = 0;
  } else {
    i = indexmap_[name];
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

Teuchos::RCP<Epetra_MultiVector> CompositeVector::ViewComponent(std::string name, bool ghosted) {
  unsigned int i;
  if (name == "" && num_components_ == 1) {
    i = 0;
  } else {
    i = indexmap_[name];
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
// Actually mutable doesn't seem to work, let's try casting away const...
Teuchos::RCP<const CompositeVector> CompositeVector::ViewOwned() const {
  if (owned_composite_ == Teuchos::null) {
    CompositeVector* temp_cv = const_cast<CompositeVector*>(this);
    return ViewOwned();
  } else {
    return owned_composite_;
  }
};

Teuchos::RCP<CompositeVector> CompositeVector::ViewOwned() {
  if (owned_composite_ == Teuchos::null) {
    // create the owned view
    owned_composite_ = Teuchos::rcp(new CompositeVector(*this, CONSTRUCT_WITHOUT_DATA));
    owned_composite_->ghosted_ = false;

    if (ghosted_) {
      for (unsigned int i = 0; i != num_components_; ++i) {
        owned_composite_->data_[i] = ViewComponent(names_[i], false);
      }
    } else {
      owned_composite_->data_ = data_;
    }
    owned_composite_->owned_data_ = owned_composite_->data_;
  }
  return owned_composite_;
};


// communicate
// -- Scatter master values to ghosted values.
// Modes shown in Epetra_CombineMode.h, but the default is Insert, which
// overwrites the current ghost value with the (unique) new master value.
void CompositeVector::ScatterMasterToGhosted(Epetra_CombineMode mode) {
#ifdef HAVE_MPI
  if (ghosted_) {
    for (unsigned int i = 0; i != num_components_; ++i) {
      // check for and create the non-ghosted view if needed
      if (owned_data_[i] == Teuchos::null) CreateOwnedView_(i);

      // check for and create the importer if needed
      if (importers_[i] == Teuchos::null) {
        if (locations_[i] == AmanziMesh::CELL) {
          const Epetra_BlockMap& target_map = mesh_->cell_map(true);
          const Epetra_BlockMap& source_map = mesh_->cell_map(false);
          importers_[i] = Teuchos::rcp(new Epetra_Import(target_map, source_map));
        } else if (locations_[i] == AmanziMesh::FACE) {
          const Epetra_BlockMap& target_map = mesh_->face_map(true);
          const Epetra_BlockMap& source_map = mesh_->face_map(false);
          importers_[i] = Teuchos::rcp(new Epetra_Import(target_map, source_map));
        } else if (locations_[i] == AmanziMesh::NODE) {
          const Epetra_BlockMap& target_map = mesh_->node_map(true);
          const Epetra_BlockMap& source_map = mesh_->node_map(false);
          importers_[i] = Teuchos::rcp(new Epetra_Import(target_map, source_map));
        } else {
          ASSERT(false);
        }
      }

      data_[i]->Import(*owned_data_[i], *importers_[i], mode);
    }
  }
#endif
};

// -- Combine ghosted values back to master values.
// Modes shown in Epetra_CombineMode.h, but the default is InsertAdd,
// where off-process values are first summed, then replace the current
// value.
void CompositeVector::GatherGhostedToMaster(Epetra_CombineMode mode) {
#ifdef HAVE_MPI
  if (ghosted_) {
    for (unsigned int i = 0; i != num_components_; ++i) {
      // check for and create the non-ghosted view if needed
      if (owned_data_[i] == Teuchos::null) CreateOwnedView_(i);

      // check for and create the importer if needed
      if (importers_[i] == Teuchos::null) {
        if (locations_[i] == AmanziMesh::CELL) {
          const Epetra_BlockMap& target_map = mesh_->cell_map(true);
          const Epetra_BlockMap& source_map = mesh_->cell_map(false);
          importers_[i] == Teuchos::rcp(new Epetra_Import(target_map, source_map));
        } else if (locations_[i] == AmanziMesh::FACE) {
          const Epetra_BlockMap& target_map = mesh_->face_map(true);
          const Epetra_BlockMap& source_map = mesh_->face_map(false);
          importers_[i] == Teuchos::rcp(new Epetra_Import(target_map, source_map));
        } else if (locations_[i] == AmanziMesh::NODE) {
          const Epetra_BlockMap& target_map = mesh_->node_map(true);
          const Epetra_BlockMap& source_map = mesh_->node_map(false);
          importers_[i] == Teuchos::rcp(new Epetra_Import(target_map, source_map));
        } else {
          ASSERT(false);
        }
      }

      owned_data_[i]->Export(*data_[i], *importers_[i], mode);
    }
  }
#endif
};

// assorted vector operations for use by time integrators?
// These may make life much easier for time integration of the full
// CompositeVector.  For instance, one could implement the BDF integrators
// doing all operations by first calling ViewOwned() to get a non-ghosted
// CompositeVector, and then using these operations.

// -- Insert value into data.
int CompositeVector::PutScalar(double scalar) {
  int ierr = 0;
  for (unsigned int i = 0; i != num_components_; ++i) {
    ierr |= data_[i]->PutScalar(scalar);
  }
  return ierr;
};

int CompositeVector::PutScalar(std::vector<double> scalar) {
  int ierr = 0;
  for (unsigned int i = 0; i != num_components_; ++i) {
    for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
      ierr = (*data_[i])(lcv_vector)->PutScalar(scalar[lcv_vector]);
      if (ierr) return ierr;
    }
  }
  return ierr;
};

int CompositeVector::PutScalar(std::string name, double scalar) {
  return data_[indexmap_[name]]->PutScalar(scalar);
};

int CompositeVector::PutScalar(std::string name, std::vector<double> scalar) {
  int ierr = 0;
  unsigned int i = indexmap_[name];
  for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
    ierr = (*data_[i])(lcv_vector)->PutScalar(scalar[lcv_vector]);
    if (ierr) return ierr;
  }
  return ierr;
};

int CompositeVector::PutScalar(std::string name, int set_id, double scalar) {
  unsigned int i = indexmap_[name];

  // valid block
  ASSERT(mesh_->valid_set_id(set_id, locations_[i]));

  // get the set
  unsigned int set_size = mesh_->get_set_size(set_id, locations_[i], AmanziMesh::OWNED);
  std::vector<unsigned int> cell_ids(set_size);
  mesh_->get_set(set_id, locations_[i], AmanziMesh::OWNED, cell_ids.begin(), cell_ids.end());

  // assign value to set
  for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
    for (std::vector<unsigned int>::iterator lcv_site = cell_ids.begin();
         lcv_site != cell_ids.end(); ++lcv_site) {
      (*(*data_[i])(lcv_vector))[*lcv_site] = scalar;
    }
  }
};

int CompositeVector::PutScalar(std::string name, int set_id, std::vector<double> scalar) {
  unsigned int i = indexmap_[name];

  // valid block
  ASSERT(mesh_->valid_set_id(set_id, locations_[i]));

  // get the set
  unsigned int set_size = mesh_->get_set_size(set_id, locations_[i], AmanziMesh::OWNED);
  std::vector<unsigned int> cell_ids(set_size);
  mesh_->get_set(set_id, locations_[i], AmanziMesh::OWNED, cell_ids.begin(), cell_ids.end());

  // assign value to set
  for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
    for (std::vector<unsigned int>::iterator lcv_site = cell_ids.begin();
         lcv_site != cell_ids.end(); ++lcv_site) {
      (*(*data_[i])(lcv_vector))[*lcv_site] = scalar[lcv_vector];
    }
  }
};

// -- this <- value*this
int CompositeVector::Scale(double value) {
  int ierr = 0;
  for (unsigned int i = 0; i != num_components_; ++i) {
    ierr = data_[i]->Scale(value);
    if (ierr) return ierr;
  }
  return ierr;
};

int CompositeVector::Scale(std::string name, double value) {
  unsigned int i = indexmap_[name];
  return data_[i]->Scale(value);
};

// -- this <- this + scalarA
int CompositeVector::Shift(double scalarA) {
  int ierr = 0;
  for (unsigned int i = 0; i != num_components_; ++i) {
    Epetra_MultiVector work(*data_[i]);
    work.PutScalar(1.0);
    ierr = data_[i]->Update(scalarA, work, 1.0);
    if (ierr) return ierr;
  }
  return ierr;
};

// -- this <- this + scalarA
int CompositeVector::Shift(std::string name, double scalarA) {
  int ierr = 0;
  unsigned int i = indexmap_[name];
  Epetra_MultiVector work(*data_[i]);
  work.PutScalar(1.0);
  ierr = data_[i]->Update(scalarA, work, 1.0);
  return ierr;
};

// -- result <- other \dot this
int CompositeVector::Dot(const CompositeVector& other, double* result) const {
  int ierr = 0;
  *result = 0.0;
  for (unsigned int i = 0; i != num_components_; ++i) {
    double intermediate_result[data_[i]->NumVectors()];
    ierr = data_[i]->Dot(*(other.data_[i]), intermediate_result);
    if (ierr) return ierr;
    for (unsigned int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
      *result += intermediate_result[lcv_vector];
    }
  }
  return ierr;
};

// -- this <- scalarA*A + scalarThis*this
CompositeVector& CompositeVector::Update(double scalarA, const CompositeVector& A, double scalarThis) {
  for (unsigned int i = 0; i != num_components_; ++i) {
    data_[i]->Update(scalarA, *A.data_[i], scalarThis);
  }
  return *this;
};

// -- this <- scalarA*A + scalarB*B + scalarThis*this
CompositeVector& CompositeVector::Update(double scalarA, const CompositeVector& A,
                 double scalarB, const CompositeVector& B, double scalarThis) {
  for (unsigned int i = 0; i != num_components_; ++i) {
    data_[i]->Update(scalarA, *A.data_[i], scalarB, *B.data_[i], scalarThis);
  }
  return *this;
};

// -- this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
int CompositeVector::Multiply(double scalarAB, const CompositeVector& A, const CompositeVector& B,
                  double scalarThis) {
  int ierr = 0;
  for (unsigned int i = 0; i != num_components_; ++i) {
    ierr = data_[i]->Multiply(scalarAB, *A.data_[i], *B.data_[i], scalarThis);
    if (ierr) return ierr;
  }
  return ierr;
};

// -- norms
int CompositeVector::NormInf(double* norm) const {
  if (norm == NULL) return 1;
  if (data_.size() == 0) return 1;

  int ierr = 0;
  *norm = 0.0;
  double norm_loc;
  for (unsigned int i = 0; i != num_components_; ++i) {
    for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
      ierr = (*data_[i])(lcv_vector)->NormInf(&norm_loc);
      if (ierr) return ierr;
      if (norm_loc > *norm) {
        *norm = norm_loc;
      }
    }
  }
  return ierr;
};

int CompositeVector::Norm1(double* norm) const {
  if (norm == NULL) return 1;
  if (data_.size() == 0) return 1;

  int ierr = 0;
  *norm = 0.0;
  double norm_loc;
  for (unsigned int i = 0; i != num_components_; ++i) {
    for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
      ierr = (*data_[i])(lcv_vector)->Norm1(&norm_loc);
      if (ierr) return ierr;
      *norm += norm_loc;
    }
  }
  return ierr;
};

int CompositeVector::Norm2(double* norm) const {
  if (norm == NULL) return 1;
  if (data_.size() == 0) return 1;

  int ierr = 0;
  *norm = 0.0;
  double norm_loc;
  for (unsigned int i = 0; i != num_components_; ++i) {
    for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
      ierr = (*data_[i])(lcv_vector)->Norm2(&norm_loc);
      if (ierr) return ierr;
      *norm += norm_loc*norm_loc;
    }
  }
  *norm = sqrt(*norm);
  return ierr;
};

void CompositeVector::Print(ostream& os) const {
  os << "Composite Vector" << std::endl;
  os << "  components: ";
  for (unsigned int i = 0; i != num_components_; ++i) {
    os << names_[i] << " ";
  }
  os << std::endl;
  for (unsigned int i = 0; i != num_components_; ++i) {
    data_[i]->Print(os);
  }
};
} // namespace

