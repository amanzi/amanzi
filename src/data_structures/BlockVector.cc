/*
  Data Structures

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Interface for BlockVector, an implementation of a slightly improved
  Epetra_MultiVector which spans multiple simplices and knows how to
  communicate itself.

  NOTE: All BlockVector data is NOT initialized to zero!
*/

#include "Epetra_MpiComm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"

#include "dbc.hh"
#include "errors.hh"

#include "BlockVector.hh"

namespace Amanzi {

// Constructor
BlockVector::BlockVector(const Comm_ptr_type& comm,
        std::vector<std::string>& names,
        std::vector<Teuchos::RCP<const Epetra_BlockMap> >& maps,
        std::vector<int> num_dofs) :
    comm_(comm),
    names_(names),
    num_dofs_(num_dofs),
    maps_(maps) {

  num_components_ = maps_.size();

  // Check consistency of input args.
  AMANZI_ASSERT(names_.size() == num_components_);
  AMANZI_ASSERT(num_dofs_.size() == num_components_);

  data_.resize(num_components_);

  // Set the size of the local portion of the map.
  sizes_.resize(num_components_);
  for (int i=0; i != num_components_; ++i) {
    indexmap_[names_[i]] = i;
    sizes_[i] = maps_[i]->NumMyElements();
  }
};


// copy constructor
BlockVector::BlockVector(const BlockVector& other) :
    comm_(other.comm_),
    indexmap_(other.indexmap_),
    num_components_(other.num_components_),
    names_(other.names_),
    num_dofs_(other.num_dofs_),
    maps_(other.maps_),
    sizes_(other.sizes_) {

  data_.resize(num_components_);
  for (int i=0; i != num_components_; ++i) {
    data_[i] = Teuchos::rcp(new Epetra_MultiVector(*other.data_[i]));
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


// Returns the length of this blockvector
int BlockVector::GlobalLength() const {
  int gl = 0;
  for (int i = 0; i != num_components_; ++i) {
    gl += data_[i]->GlobalLength();
  }
  return gl;
}


long int BlockVector::GetLocalElementCount() const {
  long int count(0);
  for (int i = 0; i != num_components_; ++i) {
    count += data_[i]->NumVectors() * data_[i]->MyLength();
  }
  return count;
}


// Check consistency of meta-data and allocate data.
void BlockVector::CreateData() {
#ifdef ENABLE_DBC
  if (data_[0] != Teuchos::null) {
    std::cout << "WARNING: CompositeVector::CreateData() called on already created vector!" << std::endl;
  }
#endif

  // create the data
  for (int i = 0; i != num_components_; ++i) {
    data_[i] = Teuchos::rcp(new Epetra_MultiVector(*maps_[i], num_dofs_[i], true));
  }
};


// View data, const version.
// I would prefer these be private, but for now...
Teuchos::RCP<const Epetra_MultiVector>
BlockVector::ViewComponent(const std::string& name) const {
  return data_[Index_(name)];
};


// View data, non-const version.
Teuchos::RCP<Epetra_MultiVector>
BlockVector::ViewComponent(const std::string& name) {
  return data_[Index_(name)];
};


// Set data
void BlockVector::SetComponent(const std::string& name,
        const Teuchos::RCP<Epetra_MultiVector>& data) {
  AMANZI_ASSERT(ComponentMap(name)->SameAs(data->Map()));
  AMANZI_ASSERT(NumVectors(name) == data->NumVectors());
  data_[Index_(name)] = data;
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
    AMANZI_ASSERT(scalar.size() == num_dofs_[i]);
    for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
      int ierr = (*data_[i])(lcv_vector)->PutScalar(scalar[lcv_vector]);
      if (ierr) return ierr;
    }
  }
  return 0;
};


// -- Insert value into component [name].
int BlockVector::PutScalar(std::string name, double scalar) {
  return data_[Index_(name)]->PutScalar(scalar);
};


// -- Insert values into data of component [name].
int BlockVector::PutScalar(std::string name, std::vector<double> scalar) {
  int i = Index_(name);
  AMANZI_ASSERT(scalar.size() == num_dofs_[i]);

  for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
    int ierr = (*data_[i])(lcv_vector)->PutScalar(scalar[lcv_vector]);
    if (ierr) return ierr;
  }
  return 0;
};


// this <- abs(this)
int BlockVector::Abs(const BlockVector& other) {
  for (int i = 0; i != num_components_; ++i) {
    int ierr = data_[i]->Abs(*other.data_[i]);
    if (ierr) return ierr;
  }
  return 0;
}


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
  return data_[Index_(name)]->Scale(value);
};


// -- this <- this + scalarA
int BlockVector::Shift(double scalarA) {
  for (int i=0; i!=num_components_; ++i) {
    Epetra_MultiVector& v = *data_[i];
    for (int j=0; j!=num_dofs_[i]; ++j) {
      for (int k=0; k!=sizes_[i]; ++k) {
        v[j][k] += scalarA;
      }
    }
  }
  return 0;
};


// Shift() applied to component name.
int BlockVector::Shift(std::string name, double scalarA) {
  int i = Index_(name);
  Epetra_MultiVector& v = *data_[i];
  for (int j=0; j!=num_dofs_[i]; ++j) {
    for (int k=0; k!=sizes_[i]; ++k) {
      v[j][k] += scalarA;
    }
  }
  return 0;
};


// this <- abs(this)
int BlockVector::Reciprocal(const BlockVector& other) {
  for (int i = 0; i != num_components_; ++i) {
    int ierr = data_[i]->Reciprocal(*other.data_[i]);
    if (ierr) return ierr;
  }
  return 0;
}


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


int BlockVector::MinValue(double* value) const {
  if (value == NULL) return 1;
  if (data_.size() == 0) return 1;

  int ierr = 0;
  double value_loc[1];

  *value = 1e+50;
  for (int i = 0; i != num_components_; ++i) {
    for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
      ierr = (*data_[i])(lcv_vector)->MinValue(value_loc);
      if (ierr) return ierr;
      *value = std::min(*value, value_loc[0]);
    }
  }
  return ierr;
};


int BlockVector::MaxValue(double* value) const {
  if (value == NULL) return 1;
  if (data_.size() == 0) return 1;

  int ierr = 0;
  double value_loc[1];

  *value = -1e+50;
  for (int i = 0; i != num_components_; ++i) {
    for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
      ierr = (*data_[i])(lcv_vector)->MaxValue(value_loc);
      if (ierr) return ierr;
      *value = std::max(*value, value_loc[0]);
    }
  }
  return ierr;
};


int BlockVector::MeanValue(double* value) const {
  if (value == NULL) return 1;
  if (data_.size() == 0) return 1;

  int ierr(0), n(0), n_loc;
  double value_loc[1];

  *value = 0.0;
  for (int i = 0; i != num_components_; ++i) {
    n_loc = data_[i]->GlobalLength(); 
    for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
      ierr = (*data_[i])(lcv_vector)->MeanValue(value_loc);
      if (ierr) return ierr;
      *value += value_loc[0] * n_loc;
      n += n_loc;
    }
  }
  *value /= n;
  return ierr;
};


// Debugging?
void BlockVector::Print(std::ostream& os, bool data_io) const {
  os << "Block Vector" << std::endl;
  os << "  components: ";
  for (int i = 0; i != num_components_; ++i) {
    os << names_[i] << "(nvec=" << data_[i]->NumVectors() << ", len=" << data_[i]->GlobalLength() << ") ";
  }
  os << std::endl;
  if (data_io) {
    for (int i = 0; i != num_components_; ++i) {
      data_[i]->Print(os);
    }
  }
};


// Populate by random numbers between -1 and 1. 
int BlockVector::Random() {
  int ierr = 0;
  for (int i = 0; i != num_components_; ++i) {
    ierr = data_[i]->Random();
    if (ierr) return ierr;
  }
  return ierr;
};


} // namespace

