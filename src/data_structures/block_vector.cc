/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for BlockVector, an implementation of a slightly improved
   Epetra_MultiVector which spans multiple simplices and knows how to
   communicate itself.

   NOTE: All BlockVector data is NOT initialized to zero!
   ------------------------------------------------------------------------- */

#include "Epetra_Vector.h"

#include "dbc.hh"
#include "errors.hh"

#include "block_vector.hh"

namespace Amanzi {

// Constructor
BlockVector::BlockVector(const Epetra_MpiComm* comm,
        std::vector<std::string>& names,
        std::vector<Teuchos::RCP<const Epetra_Map> >& maps,
        std::vector<int> num_dofs) :

    comm_(comm),
    names_(names),
    maps_(maps),
    num_dofs_(num_dofs) {

  num_components_ = maps_.size();

  // Check consistency of input args.
  ASSERT(names_.size() == num_components_);
  ASSERT(num_dofs_.size() == num_components_);

  data_.resize(num_components_);

  // Set the size of the local portion of the map.
  sizes_.resize(num_components_);
  for (int i=0; i != num_components_; ++i) {
    indexmap_[names_[i]] = i;
    sizes_[i] = maps_[i]->NumMyElements();
  }
};


// copy constructor
BlockVector::BlockVector(const BlockVector& other, ConstructMode mode) :
    comm_(other.comm_),
    names_(other.names_),
    maps_(other.maps_),
    num_dofs_(other.num_dofs_),
    num_components_(other.num_components_),
    sizes_(other.sizes_),
    indexmap_(other.indexmap_) {

  data_.resize(num_components_);

  if (mode == CONSTRUCT_WITH_NEW_DATA) {
    for (int i=0; i != num_components_; ++i) {
      data_[i] = Teuchos::rcp(new Epetra_MultiVector(*other.data_[i]));
    }
  } else if (mode == CONSTRUCT_WITH_OLD_DATA) {
    for (int i=0; i != num_components_; ++i) {
      data_[i] = other.data_[i];
    }
  }
};


// Assigment.
BlockVector& BlockVector::operator=(const BlockVector& other) {
  if (this != &other) {

#ifdef ENABLE_DBC
    if (num_components_ != other.num_components_) {
      Errors::Message message("Attempted assignment of non-compatible BlockVectors.");
      Exceptions::amanzi_throw(message);
    }

    for (int i = 0; i != num_components_; ++i) {
      if ((names_[i] != other.names_[i]) ||
          (num_dofs_[i] != other.num_dofs_[i]) ||
          (maps_[i] != other.maps_[i])) {
        Errors::Message
          message("Attempted assignment of non-compatible BlockVectors.");
        Exceptions::amanzi_throw(message);
      }
    }
#endif // ENABLE_DBC

    for (int i = 0; i != num_components_; ++i) {
      if (other.data_[i] == Teuchos::null) {
        data_[i] = Teuchos::null;
      } else if (data_[i] == Teuchos::null) {
        data_[i] = Teuchos::rcp(new Epetra_MultiVector(*other.data_[i]));
      } else {
        *data_[i] = *other.data_[i];
      }
    }
  }
  return *this;
};


// Check consistency of meta-data and allocate data.
void BlockVector::CreateData() {
  // create the data
  for (int i = 0; i != num_components_; ++i) {
    data_[i] = Teuchos::rcp(new Epetra_MultiVector(*maps_[i], num_dofs_[i], false));
  }
};


// View data, const version.
// I would prefer these be private, but for now...
Teuchos::RCP<const Epetra_MultiVector>
BlockVector::ViewComponent(std::string name) const {
  return data_[index_(name)];
};


// View data, non-const version.
Teuchos::RCP<Epetra_MultiVector>
BlockVector::ViewComponent(std::string name) {
  return data_[index_(name)];
};


// Set data
void BlockVector::SetComponent(std::string name,
        const Teuchos::RCP<Epetra_MultiVector>& data) {
  ASSERT(map(name)->SameAs(data->Map()));
  ASSERT(num_dofs(name) == data->NumVectors());
  data_[index_(name)] = data;
}

// Vector operations.
// -- Insert value into data.
int BlockVector::PutScalar(double scalar) {
  int ierr = 0;
  for (int i = 0; i != num_components_; ++i) {
    ierr = data_[i]->PutScalar(scalar);
    if (ierr) return ierr;
  }
  return ierr;
};


// -- Insert values into data, by DOF, not by component!
int BlockVector::PutScalar(std::vector<double> scalar) {
  for (int i = 0; i != num_components_; ++i) {
    ASSERT(scalar.size() == num_dofs_[i]);
    for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
      int ierr = (*data_[i])(lcv_vector)->PutScalar(scalar[lcv_vector]);
      if (ierr) return ierr;
    }
  }
  return 0;
};


// -- Insert value into component [name].
int BlockVector::PutScalar(std::string name, double scalar) {
  return data_[index_(name)]->PutScalar(scalar);
};


// -- Insert values into data of component [name].
int BlockVector::PutScalar(std::string name, std::vector<double> scalar) {
  int i = index_(name);
  ASSERT(scalar.size() == num_dofs_[i]);

  for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
    int ierr = (*data_[i])(lcv_vector)->PutScalar(scalar[lcv_vector]);
    if (ierr) return ierr;
  }
  return 0;
};


// -- this <- value*this
int BlockVector::Scale(double value) {
  for (int i = 0; i != num_components_; ++i) {
    int ierr = data_[i]->Scale(value);
    if (ierr) return ierr;
  }
  return 0;
};


// Scale() applied to component name.
int BlockVector::Scale(std::string name, double value) {
  return data_[index_(name)]->Scale(value);
};


// -- this <- this + scalarA
int BlockVector::Shift(double scalarA) {
  for (int i=0; i!=num_components_; ++i) {
    for (int j=0; j!=num_dofs_[i]; ++j) {
      for (int k=0; k!=sizes_[i]; ++k) {
        operator()(names_[i], j, k) += scalarA;
      }
    }
  }
  return 0;
};


// Shift() applied to component name.
int BlockVector::Shift(std::string name, double scalarA) {
  int i = index_(name);
  for (int j=0; j!=num_dofs_[i]; ++j) {
    for (int k=0; k!=sizes_[i]; ++k) {
      operator()(names_[i], j, k) += scalarA;
    }
  }
  return 0;
};


// -- result <- other \dot this
int BlockVector::Dot(const BlockVector& other, double* result) const {
  *result = 0.0;
  for (int i = 0; i != num_components_; ++i) {
    double intermediate_result[data_[i]->NumVectors()];
    int ierr = data_[i]->Dot(*(other.data_[i]), intermediate_result);
    if (ierr) return ierr;
    for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
      *result += intermediate_result[lcv_vector];
    }
  }
  return 0;
};


// -- this <- scalarA*A + scalarThis*this
BlockVector& BlockVector::Update(double scalarA, const BlockVector& A, double scalarThis) {
  for (int i = 0; i != num_components_; ++i) {
    data_[i]->Update(scalarA, *A.data_[i], scalarThis);
  }
  return *this;
};


// -- this <- scalarA*A + scalarB*B + scalarThis*this
BlockVector& BlockVector::Update(double scalarA, const BlockVector& A,
                 double scalarB, const BlockVector& B, double scalarThis) {
  for (int i = 0; i != num_components_; ++i) {
    data_[i]->Update(scalarA, *A.data_[i], scalarB, *B.data_[i], scalarThis);
  }
  return *this;
};


// -- this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
int BlockVector::Multiply(double scalarAB, const BlockVector& A, const BlockVector& B,
                  double scalarThis) {
  int ierr = 0;
  for (int i = 0; i != num_components_; ++i) {
    ierr = data_[i]->Multiply(scalarAB, *A.data_[i], *B.data_[i], scalarThis);
    if (ierr) return ierr;
  }
  return ierr;
};

// -- this <- scalarAB * B / A + scalarThis*this  (/ is the elementwise division
int BlockVector::ReciprocalMultiply(double scalarAB, const BlockVector& A, const BlockVector& B,
                  double scalarThis) {
  int ierr = 0;
  for (int i = 0; i != num_components_; ++i) {
    ierr = data_[i]->ReciprocalMultiply(scalarAB, *A.data_[i], *B.data_[i], scalarThis);
    if (ierr) return ierr;
  }
  return ierr;
};


// -- norms
int BlockVector::NormInf(double* norm) const {
  if (norm == NULL) return 1;
  if (data_.size() == 0) return 1;

  int ierr = 0;
  *norm = 0.0;
  double norm_loc;
  for (int i = 0; i != num_components_; ++i) {
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


int BlockVector::Norm1(double* norm) const {
  if (norm == NULL) return 1;
  if (data_.size() == 0) return 1;

  int ierr = 0;
  *norm = 0.0;
  double norm_loc;
  for (int i = 0; i != num_components_; ++i) {
    for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
      ierr = (*data_[i])(lcv_vector)->Norm1(&norm_loc);
      if (ierr) return ierr;
      *norm += norm_loc;
    }
  }
  return ierr;
};


int BlockVector::Norm2(double* norm) const {
  if (norm == NULL) return 1;
  if (data_.size() == 0) return 1;

  int ierr = 0;
  *norm = 0.0;
  double norm_loc;
  for (int i = 0; i != num_components_; ++i) {
    for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
      ierr = (*data_[i])(lcv_vector)->Norm2(&norm_loc);
      if (ierr) return ierr;
      *norm += norm_loc*norm_loc;
    }
  }
  *norm = sqrt(*norm);
  return ierr;
};


// Debugging?
void BlockVector::Print(ostream& os) const {
  os << "Block Vector" << std::endl;
  os << "  components: ";
  for (int i = 0; i != num_components_; ++i) {
    os << names_[i] << " ";
  }
  os << std::endl;
  for (int i = 0; i != num_components_; ++i) {
    data_[i]->Print(os);
  }
};

} // namespace

