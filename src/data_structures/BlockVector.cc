/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

/*

  BlockVector

  Interface for BlockVector, an implementation of a slightly improved
  Epetra_MultiVector which spans multiple simplices and knows how to
  communicate itself.

*/
  
#include <numeric>
#include "AmanziComm.hh"
#include "AmanziMap.hh"
#include "AmanziVector.hh"

#include "dbc.hh"
#include "errors.hh"

#include "BlockVector.hh"

namespace Amanzi {

// Constructor
BlockVector::BlockVector(const Comm_ptr_type& comm,
        const std::vector<std::string>& names,
        const std::vector<Map_ptr_type >& maps,
        const std::vector<int> num_dofs) :
    names_(names),
    maps_(maps),
    num_dofs_(num_dofs),
    comm_(comm) {

  num_components_ = maps_.size();

  // Check consistency of input args.
  AMANZI_ASSERT(names_.size() == num_components_);
  AMANZI_ASSERT(num_dofs_.size() == num_components_);

  data_.resize(num_components_);

  // Set the size of the local portion of the map.
  sizes_.resize(num_components_);
  for (int i=0; i != num_components_; ++i) {
    indexmap_[names_[i]] = i;
    sizes_[i] = maps_[i]->getNodeNumElements();
  }
};


// copy constructor
BlockVector::BlockVector(const BlockVector& other) :
    comm_(other.comm_),
    names_(other.names_),
    maps_(other.maps_),
    num_dofs_(other.num_dofs_),
    num_components_(other.num_components_),
    sizes_(other.sizes_),
    indexmap_(other.indexmap_) {

  data_.resize(num_components_);
  for (int i=0; i != num_components_; ++i) {
    data_[i] = Teuchos::rcp(new MultiVector_type(*other.data_[i], Teuchos::Copy));
  }
};


// // Assigment.
// BlockVector& BlockVector::operator=(const BlockVector& other) {
//   if (this != &other) {

// #ifdef ENABLE_DBC
//     if (num_components_ != other.num_components_) {
//       Errors::Message message("Attempted assignment of non-compatible BlockVectors.");
//       Exceptions::amanzi_throw(message);
//     }

//     for (int i = 0; i != num_components_; ++i) {
//       if ((names_[i] != other.names_[i]) ||
//           (num_dofs_[i] != other.num_dofs_[i]) ||
//           (maps_[i] != other.maps_[i])) {
//         Errors::Message message("Attempted assignment of non-compatible BlockVectors.");
//         Exceptions::amanzi_throw(message);
//       }
//     }
// #endif // ENABLE_DBC

//     for (int i = 0; i != num_components_; ++i) {
//       if (other.data_[i] == Teuchos::null) {
//         data_[i] = Teuchos::null;
//       } else if (data_[i] == Teuchos::null) {
//         data_[i] = Teuchos::rcp(new MultiVector_type(*other.data_[i], Teuchos::Copy));
//       } else {
//         *data_[i] = *other.data_[i];
//       }
//     }
//   }
//   return *this;
// };


// Returns the length of this blockvector
int BlockVector::GlobalLength() const {
  int gl = 0;
  for (int i = 0; i != num_components_; ++i) {
    gl += data_[i]->getGlobalLength();
  }
  return gl;
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
    data_[i] = Teuchos::rcp(new MultiVector_type(maps_[i], num_dofs_[i], true));
  }
};


// View data, const version.
// I would prefer these be private, but for now...
Teuchos::RCP<const MultiVector_type>
BlockVector::GetComponent(std::string name) const {
  return data_[Index_(name)];
}


// View data, non-const version.
MultiVector_ptr_type
BlockVector::GetComponent(std::string name) {
  return data_[Index_(name)];
};


// Set data
void BlockVector::SetComponent(std::string name,
                                 const MultiVector_ptr_type& data) {
  AMANZI_ASSERT(ComponentMap(name)->isSameAs(*data->getMap()));
  AMANZI_ASSERT(NumVectors(name) == data->getNumVectors());
  data_[Index_(name)] = data;
}

// Vector operations.
// -- Insert value into data.
int BlockVector::PutScalar(double scalar) {
  int ierr = 0;
  for (int i = 0; i != num_components_; ++i) {
    data_[i]->putScalar(scalar);
  }
  return ierr;
};


// // -- Insert values into data, by DOF, not by component!
// int BlockVector::PutScalar(std::vector<double> scalar) {
//   for (int i = 0; i != num_components_; ++i) {
//     AMANZI_ASSERT(scalar.size() == num_dofs_[i]);
//     for (int lcv_vector = 0; lcv_vector != data_[i]->getNumVectors(); ++lcv_vector) {
//       (*data_[i])(lcv_vector)->putScalar(scalar[lcv_vector]);
//     }
//   }
//   return 0;
// };


// -- Insert value into component [name].
int BlockVector::PutScalar(std::string name, double scalar) {
  data_[Index_(name)]->putScalar(scalar);
  return 0;
};


// // -- Insert values into data of component [name].
// int BlockVector::PutScalar(std::string name, std::vector<double> scalar) {
//   int i = Index_(name);
//   AMANZI_ASSERT(scalar.size() == num_dofs_[i]);

//   for (int lcv_vector = 0; lcv_vector != data_[i]->getNumVectors(); ++lcv_vector) {
//     (*data_[i])(lcv_vector)->putScalar(scalar[lcv_vector]);
//   }
//   return 0;
// };


// this <- abs(this)
int BlockVector::Abs(const BlockVector& other) {
  for (int i = 0; i != num_components_; ++i) {
    data_[i]->abs(*other.data_[i]);
  }
  return 0;
}


// -- this <- value*this
int BlockVector::Scale(double value) {
  for (int i = 0; i != num_components_; ++i) {
    data_[i]->scale(value);
  }
  return 0;
};


// Scale() applied to component name.
int BlockVector::Scale(std::string name, double value) {
  data_[Index_(name)]->scale(value);
  return 0;
};


// this <- abs(this)
int BlockVector::Reciprocal(const BlockVector& other) {
  for (int i = 0; i != num_components_; ++i) {
    data_[i]->reciprocal(*other.data_[i]);
  }
  return 0;
}


// -- result <- other \dot this
int BlockVector::Dot(const BlockVector& other, double* result) const {
  *result = 0.0;
  for (int i = 0; i != num_components_; ++i) {
    Teuchos::Array<double> intermediate_result(data_[i]->getNumVectors());
    data_[i]->dot(*(other.data_[i]), intermediate_result());
    for (int lcv_vector = 0; lcv_vector != data_[i]->getNumVectors(); ++lcv_vector) {
      *result += intermediate_result[lcv_vector];
    }
  }
  return 0;
};


// -- this <- scalarA*A + scalarThis*this
BlockVector& BlockVector::Update(double scalarA, const BlockVector& A, double scalarThis) {
  for (int i = 0; i != num_components_; ++i) {
    data_[i]->update(scalarA, *A.data_[i], scalarThis);
  }
  return *this;
};


// -- this <- scalarA*A + scalarB*B + scalarThis*this
BlockVector& BlockVector::Update(double scalarA, const BlockVector& A,
                 double scalarB, const BlockVector& B, double scalarThis) {
  for (int i = 0; i != num_components_; ++i) {
    data_[i]->update(scalarA, *A.data_[i], scalarB, *B.data_[i], scalarThis);
  }
  return *this;
};


// -- this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
int BlockVector::Multiply(double scalarAB, const BlockVector& A, const BlockVector& B,
                  double scalarThis) {
  for (int i = 0; i != num_components_; ++i) {
    if (A.data_[i]->getNumVectors() > 1) {
      Errors::Message message("Not implemented multiply: Tpetra does not provide elementwise multiply.");
      Exceptions::amanzi_throw(message);
    }
    data_[i]->elementWiseMultiply(scalarAB, *A.data_[i]->getVector(0), *B.data_[i], scalarThis);
  }
  return 0;
};

// // -- this <- scalarAB * B / A + scalarThis*this  (/ is the elementwise division
// int BlockVector::ReciprocalMultiply(double scalarAB, const BlockVector& A, const BlockVector& B,
//                   double scalarThis) {
//   for (int i = 0; i != num_components_; ++i) {
//     data_[i]->reciprocalMultiply(scalarAB, *A.data_[i], *B.data_[i], scalarThis);
//   }
//   return 0;
// };


// -- norms
int BlockVector::NormInf(double* norm) const {
  if (norm == NULL) return 1;
  if (data_.size() == 0) return 1;

  *norm = 0.0;
  for (int i = 0; i != num_components_; ++i) {
    Teuchos::Array<double> norm_locs(data_[i]->getNumVectors());
    data_[i]->normInf(norm_locs());
    double my_norm_loc = *std::max_element(norm_locs.begin(), norm_locs.end());
    *norm = std::max(my_norm_loc, *norm);
  }
  return 0;
};


int BlockVector::Norm1(double* norm) const {
  if (norm == NULL) return 1;
  if (data_.size() == 0) return 1;

  *norm = 0.0;
  for (int i = 0; i != num_components_; ++i) {
    Teuchos::Array<double> norm_locs(data_[i]->getNumVectors());
    data_[i]->norm1(norm_locs());
    *norm += std::accumulate(norm_locs.begin(), norm_locs.end(), 0);
  }
  return 0;
};


int BlockVector::Norm2(double* norm) const {
  if (norm == NULL) return 1;
  if (data_.size() == 0) return 1;

  *norm = 0.0;
  for (int i = 0; i != num_components_; ++i) {
    Teuchos::Array<double> norm_locs(data_[i]->getNumVectors());
    data_[i]->norm2(norm_locs());
    *norm += std::accumulate(norm_locs.begin(), norm_locs.end(), 0, [](double x, double y) { return x + y*y; } );
  }
  *norm = sqrt(*norm);
  return 0;
};



// int BlockVector::MinValue(double* value) const {
//   if (value == NULL) return 1;
//   if (data_.size() == 0) return 1;

//   int ierr = 0;
//   double value_loc[1];

//   *value = 1e+50;
//   for (int i = 0; i != num_components_; ++i) {
//     for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
//       ierr = (*data_[i])(lcv_vector)->MinValue(value_loc);
//       if (ierr) return ierr;
//       *value = std::min(*value, value_loc[0]);
//     }
//   }
//   return ierr;
// };


// int BlockVector::MaxValue(double* value) const {
//   if (value == NULL) return 1;
//   if (data_.size() == 0) return 1;

//   int ierr = 0;
//   double value_loc[1];

//   *value = -1e+50;
//   for (int i = 0; i != num_components_; ++i) {
//     for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
//       ierr = (*data_[i])(lcv_vector)->MaxValue(value_loc);
//       if (ierr) return ierr;
//       *value = std::max(*value, value_loc[0]);
//     }
//   }
//   return ierr;
// };


// int BlockVector::MeanValue(double* value) const {
//   if (value == NULL) return 1;
//   if (data_.size() == 0) return 1;

//   int ierr(0), n(0), n_loc;
//   double value_loc[1];

//   *value = 0.0;
//   for (int i = 0; i != num_components_; ++i) {
//     n_loc = data_[i]->GlobalLength(); 
//     for (int lcv_vector = 0; lcv_vector != data_[i]->NumVectors(); ++lcv_vector) {
//       ierr = (*data_[i])(lcv_vector)->MeanValue(value_loc);
//       if (ierr) return ierr;
//       *value += value_loc[0] * n_loc;
//       n += n_loc;
//     }
//   }
//   *value /= n;
//   return ierr;
// };


// Debugging?
void BlockVector::Print(std::ostream& os, bool data_io) const {
  os << "Block Vector" << std::endl;
  os << "  components: ";
  for (int i = 0; i != num_components_; ++i) {
    os << names_[i] << "(" << data_[i]->getNumVectors() << ") ";
  }
  os << std::endl;
  if (data_io) {
    for (int i = 0; i != num_components_; ++i) {
      data_[i]->print(os);
    }
  }
};


// Populate by random numbers between -1 and 1. 
int BlockVector::Random() {
  for (int i = 0; i != num_components_; ++i) {
    data_[i]->randomize();
  }
  return 0;
};


} // namespace

