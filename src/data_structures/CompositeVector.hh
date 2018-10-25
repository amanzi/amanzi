/* -------------------------------------------------------------------------
ATS & Amanzi

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon (ecoon@lanl.gov)

Interface for CompositeVector, an implementation of a block Tpetra
vector which spans multiple simplices and knows how to communicate
itself.

CompositeVectors are a collection of vectors defined on a common mesh and
communicator.  Each vector, or component, has a name (used as a key), a mesh
Entity_kind (CELL, FACE, NODE, or BOUNDARY_FACE), and a number of degrees of
freedom (dofs).  This, along with the Map_type provided from the mesh on a
given Entity_kind, is enough to create a Vector

Note that construction of the CompositeVector does not allocate the
Tpetra_Vector.  CreateData() must be called before usage.

Access using operator() is slow, and should only be used for debugging.
Prefer to use the ViewComponent() accessors.

This vector provides the duck-type interface Vec and may be used with time
integrators/nonlinear solvers.

DOCUMENT VANDELAY HERE! FIX ME --etc
------------------------------------------------------------------------- */

#ifndef AMANZI_COMPOSITEVECTOR_HH_
#define AMANZI_COMPOSITEVECTOR_HH_

#define CV_ENABLE_SET_FROM_OPERATOR 0

#include <vector>
#include "Teuchos_RCP.hpp"
#include "AmanziTypes.hh"

#include "dbc.hh"
#include "Mesh.hh"
#include "data_structures_types.hh"
#include "BlockVector.hh"
#include "CompositeVectorSpace.hh"

namespace Amanzi {

class CompositeVector {

public:
  // -- Constructors --

  // Constructor from a CompositeVectorSpace (which is like a VectorSpace).
  CompositeVector(const CompositeVectorSpace& space);
  CompositeVector(const CompositeVectorSpace& space, bool ghosted);

  // Copy constructor.
  CompositeVector(const CompositeVector& other, InitMode mode=INIT_MODE_COPY);
  CompositeVector(const CompositeVector& other, bool ghosted,
                  InitMode mode=INIT_MODE_COPY);

  // Assignment operator.
  CompositeVector& operator=(const CompositeVector& other);

  // -- Accessors to meta-data --

  // Space/VectorSpace/Map accessor.
  const CompositeVectorSpace& Map() const { return *map_; }

  // CompositeVector maintains its own ghosted value.
  bool Ghosted() const { return ghosted_; }

  // Much accessor functionality is delegated to the VectorSpace.
  typedef std::vector<std::string>::const_iterator name_iterator;
  name_iterator begin() const { return map_->begin(); }
  name_iterator end() const { return map_->end(); }
  unsigned int size() const { return map_->size(); }
  Comm_ptr_type Comm() const { return map_->Comm(); }
  Teuchos::RCP<const AmanziMesh::Mesh> Mesh() const { return map_->Mesh(); }
  bool HasComponent(std::string name) const { return map_->HasComponent(name); }
  int NumComponents() const { return size(); }
  int NumVectors(std::string name) const { return map_->NumVectors(name); }
  int GlobalLength() { return mastervec_->GlobalLength(); }
  AmanziMesh::Entity_kind Location(std::string name) const { return map_->Location(name); }

  // Provides the size of each component's vector, either ghosted or non-ghosted.
  unsigned int size(std::string name, bool ghosted=false) const {
    return ghosted ? ghostvec_->size(name) : mastervec_->size(name); }

  int GlobalLength() const { return mastervec_->GlobalLength(); }
  
  // Access the VectorSpace for each component.
  Teuchos::RCP<const Map_type> ComponentMap(std::string name,
          bool ghosted=false) const {
    return ghosted ? ghostvec_->ComponentMap(name) : mastervec_->ComponentMap(name);
  }

  // -- View data. --
  // Access a view of a single component's data.
  //
  // Const access -- this does not tag as changed.
  Teuchos::RCP<const Epetra_MultiVector>
  ViewComponent(std::string name, bool ghosted=false) const;

  // View entries in the vectors
  //
  // Return-by-value, this does not tag as changed.
  double operator()(std::string name, int i, int j) const {
    return (*ghostvec_)(name,i,j);
  }
  double operator()(std::string name, int j) const {
    return (*ghostvec_)(name,0,j);
  }

  // -- Set data. --

  // Access a view of a single component's data.
  //
  // Non-const access -- tags changed.
  Teuchos::RCP<Epetra_MultiVector>
  ViewComponent(std::string name, bool ghosted=false);

#if CV_ENABLE_SET_FROM_OPERATOR
  // Set entries in the vectors.
  //
  // Using these is VERY STRONGLY DISCOURAGED.  Instead, call ViewComponent()
  // and set entries in the MultiVector.  THESE ARE VERY VERY SLOW (But they
  // can be handy in debugging.)  Tags changed.
  double& operator()(std::string name, int i, int j) {
    return (*ghostvec_)(name,i,j);
  }

  // Set entries in the 0th vector.
  //
  // Using these is VERY STRONGLY DISCOURAGED.  Instead, call ViewComponent()
  // and set entries in the MultiVector.  THESE ARE VERY VERY SLOW (But they
  // can be handy in debugging.)  Tags changed.
  double& operator()(std::string name, int j) {
    return (*ghostvec_)(name,0,j);
  }
#endif


  // Set block by pointer if possible, copy if not?????? FIX ME --etc
  void SetComponent(std::string name, const Teuchos::RCP<Epetra_MultiVector>& data);

  // Scatter master values to ghosted values, on all components (INSERT mode).
  //
  // Note that although scatter changes things, it doesn't change master
  // data, so we allow it to work on const.  This is necessary for a
  // non-owning PK to communicate a non-owned vector.
  //
  // Note that the scatter is ifneeded, unless force=true.
  void ScatterMasterToGhosted(bool force=false) const;

  // Scatter master values to ghosted values, on one components (INSERT mode).
  //
  // Modes shown in Epetra_CombineMode.h, but the default is Insert, which
  // overwrites the current ghost value with the (unique) new master value.
  //
  // Note that although scatter changes things, it doesn't change master
  // data, so we allow it to work on const.  This is necessary for a
  // non-owning PK to communicate a non-owned vector.
  //
  // Note that the scatter is ifneeded, unless force=true.
  void ScatterMasterToGhosted(std::string name, bool force=false) const;
  void ScatterMasterToGhosted(const char* name, bool force=false) const {
    ScatterMasterToGhosted(std::string(name), force); }

  // Scatter master values to ghosted values, on all components, in a mode.
  //
  // Modes shown in Epetra_CombineMode.h, but the default is Insert, which
  // overwrites the current ghost value with the (unique) new master value.
  //
  // Note that although scatter changes things, it doesn't change master
  // data, so we allow it to work on const.  This is necessary for a
  // non-owning PK to communicate a non-owned vector.
  //
  // This Scatter() is not managed, and is always done.  Tags changed.
  void ScatterMasterToGhosted(Epetra_CombineMode mode) const;

  // Scatter master values to ghosted values, on all components, in a mode.
  //
  // Modes shown in Epetra_CombineMode.h, but the default is Insert, which
  // overwrites the current ghost value with the (unique) new master value.
  //
  // Note that although scatter changes things, it doesn't change master
  // data, so we allow it to work on const.  This is necessary for a
  // non-owning PK to communicate a non-owned vector.
  //
  // This Scatter() is not managed, and is always done.  Tags changed.
  void ScatterMasterToGhosted(std::string name, Epetra_CombineMode mode) const;

  // Combine ghosted values back to master values.
  //
  // Modes shown in Epetra_CombineMode.h, but the default is Add,
  // where off-process values are first summed into the on-process value.
  //
  // This Scatter() is not managed, and is always done.  Tags changed.
  void GatherGhostedToMaster(Epetra_CombineMode mode=Add);
  void GatherGhostedToMaster(std::string name, Epetra_CombineMode mode=Add);

  // returns non-empty importer
  const Teuchos::RCP<Epetra_Import>& importer(std::string name);

  // -- Assorted vector operations, this implements a Vec --

  // Sets all vectors to value.
  int PutScalar(double scalar);

  // Sets all vectors to value including ghosted elements.
  // Different name is given so it cannot be used in a templated code.   
  int PutScalarMasterAndGhosted(double scalar);

  // Sets ghost elements to value.
  // Different name is given so it cannot be used in a templated code.   
  int PutScalarGhosted(double scalar);

  // v(name,:,:) = scalar
  int PutScalar(std::string name, double scalar);

  // v(name,i,:) = scalar[i]
  int PutScalar(std::string name, std::vector<double> scalar);

  // this <- scalar*this
  int Scale(double scalar);
  int ScaleMasterAndGhosted(double scalar);

  // this <- abs(this)
  int Abs(const CompositeVector& other);
  
  // this(name,:,:) <- scalar*this(name,:,:)
  int Scale(std::string name, double scalar);

  // this <- this + scalar
  int Shift(double scalar);

  // this(name,:,:) <- scalar + this(name,:,:)
  int Shift(std::string name, double scalar);

  // this <- element wise reciprocal(this)
  int Reciprocal(const CompositeVector& other);
  
  // result <- other \dot this
  int Dot(const CompositeVector& other, double* result) const;

  // this <- scalarA*A + scalarThis*this
  CompositeVector& Update(double scalarA, const CompositeVector& A, double scalarThis);

  // this <- scalarA*A + scalarB*B + scalarThis*this
  CompositeVector& Update(double scalarA, const CompositeVector& A,
                          double scalarB, const CompositeVector& B, double scalarThis);

  // this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
  int Multiply(double scalarAB, const CompositeVector& A, const CompositeVector& B,
               double scalarThis);

  // this <- scalarAB * A^-1@B + scalarThis*this  (@ is the elementwise product
  int ReciprocalMultiply(double scalarAB, const CompositeVector& A,
                         const CompositeVector& B, double scalarThis);

  // -- norms --
  int NormInf(double* norm) const;
  int Norm1(double* norm) const;
  int Norm2(double* norm) const;

  int MinValue(double* value) const;
  int MaxValue(double* value) const;
  int MeanValue(double* value) const;

  void MinValue(std::map<std::string, double>& value) const;
  void MaxValue(std::map<std::string, double>& value) const;
  void MeanValue(std::map<std::string, double>& value) const;

  // -- Utilities --

  // Write components to outstream.
  void Print(std::ostream &os, bool data_io = true) const;

  // Populate by random numbers between -1 and 1.
  int Random();

 private:
  void InitMap_(const CompositeVectorSpace& space);
  void InitData_(const CompositeVector& other, InitMode mode);
  void CreateData_();

  int Index_(std::string name) const {
    std::map<std::string, int>::const_iterator item = indexmap_.find(name);
    AMANZI_ASSERT(item != indexmap_.end());
    return item->second;
  }


  // The Vandelay is an Importer/Exporter which allows face unknowns
  // to be spoofed as boundary face unknowns.
  void CreateVandelay_() const;
  void ApplyVandelay_() const;


 private:
  Teuchos::RCP<const CompositeVectorSpace> map_;
  bool ghosted_;

  // data enumerating the blocks
  std::map< std::string, int > indexmap_;
  std::vector<std::string> names_;

  // data containing the blocks
  mutable Teuchos::RCP<BlockVector> ghostvec_;
  Teuchos::RCP<BlockVector> mastervec_;

  // importers for scatter/gather operation
  mutable std::vector<Teuchos::RCP<Epetra_Import> > importers_;

  // importer and vector for boundary data
  mutable Teuchos::RCP<Epetra_Import> vandelay_importer_;
  mutable Teuchos::RCP<Epetra_MultiVector> vandelay_vector_;
};


inline int
CompositeVector::PutScalar(double scalar) {
  return mastervec_->PutScalar(scalar);
}

inline int
CompositeVector::PutScalarMasterAndGhosted(double scalar) {
  return ghostvec_->PutScalar(scalar);
}

inline int
CompositeVector::PutScalarGhosted(double scalar) {
  for (int lcv_comp = 0; lcv_comp != NumComponents(); ++lcv_comp) {
    int size_owned = mastervec_->size(names_[lcv_comp]);
    int size_ghosted = ghostvec_->size(names_[lcv_comp]);

    Epetra_MultiVector& vec = *ghostvec_->ViewComponent(names_[lcv_comp]);
    for (int j = 0; j != vec.NumVectors(); ++j) {
      for (int i = size_owned; i != size_ghosted; ++i) {
        vec[j][i] = scalar;
      }
    }
  }
  return 0;
}

inline int
CompositeVector::PutScalar(std::string name, double scalar) {
  return mastervec_->PutScalar(name, scalar);
}

inline int
CompositeVector::PutScalar(std::string name, std::vector<double> scalar) {
  return mastervec_->PutScalar(name, scalar);
}

inline int
CompositeVector::Abs(const CompositeVector& other) {
  return mastervec_->Abs(*other.mastervec_);
}

inline int
CompositeVector::Scale(double scalar) {
  return mastervec_->Scale(scalar);
}

inline int
CompositeVector::ScaleMasterAndGhosted(double scalar) {
  return ghostvec_->Scale(scalar);
}

inline int
CompositeVector::Scale(std::string name, double scalar) {
  return mastervec_->Scale(name, scalar);
}

inline int
CompositeVector::Shift(double scalar) {
  return mastervec_->Shift(scalar);
}

inline int
CompositeVector::Shift(std::string name, double scalar) {
  return mastervec_->Shift(name, scalar);
}

inline int
CompositeVector::Reciprocal(const CompositeVector& other) {
  return mastervec_->Reciprocal(*other.mastervec_);
}

inline int
CompositeVector::NormInf(double* norm) const {
  return mastervec_->NormInf(norm);
}

inline int
CompositeVector::Norm1(double* norm) const {
  return mastervec_->Norm1(norm);
}

inline int
CompositeVector::Norm2(double* norm) const {
  return mastervec_->Norm2(norm);
}

inline int
CompositeVector::MinValue(double* value) const {
  return mastervec_->MinValue(value);
}

inline int
CompositeVector::MaxValue(double* value) const {
  return mastervec_->MaxValue(value);
}

inline int
CompositeVector::MeanValue(double* value) const {
  return mastervec_->MeanValue(value);
}

inline void
CompositeVector::Print(std::ostream &os, bool data_io) const {
  return mastervec_->Print(os, data_io);
}

inline int
CompositeVector::Random() {
  return mastervec_->Random();
}


// -----------------------------------------------------------------------------
// Non-member functions.
// -----------------------------------------------------------------------------
void DeriveFaceValuesFromCellValues(CompositeVector&);

void AddComponent(Teuchos::RCP<CompositeVector> cv,
                  const std::string& name, AmanziMesh::Entity_kind kind, int dim);

} // namespace

#endif
