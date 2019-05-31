/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

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

*/

#ifndef AMANZI_COMPOSITEVECTOR_DECL_HH_
#define AMANZI_COMPOSITEVECTOR_DECL_HH_

#include <vector>
#include "Teuchos_RCP.hpp"
#include "AmanziTypes.hh"
#include "AmanziVector.hh"

#include "dbc.hh"
#include "CompVector.hh"
#include "CompositeMap.hh"

namespace Amanzi {

//
// Class interface
//
template<typename Scalar>
class CompositeVector_ {

public:
  // -- Constructors --
  // Constructor from a CompositeMap (which is like a Map).
  CompositeVector_(const Teuchos::RCP<const CompositeMap>& space);
  CompositeVector_(const Teuchos::RCP<const CompositeMap>& space, bool ghosted);

  // Copy constructor.
  CompositeVector_(const CompositeVector_<Scalar>& other, InitMode mode=INIT_MODE_COPY);
  CompositeVector_(const CompositeVector_<Scalar>& other, bool ghosted,
                  InitMode mode=INIT_MODE_COPY);

  ~CompositeVector_() = default;

  // Assignment operator.
  CompositeVector_<Scalar>& operator=(const CompositeVector_<Scalar>& other);
  
  // -- Accessors to meta-data --
  // Space/VectorSpace/Map accessor.
  Teuchos::RCP<const CompositeMap> Map() const { return map_; }

  // Much accessor functionality is delegated to the VectorSpace.
  using name_iterator = CompositeMap::name_iterator;
  name_iterator begin() const;
  name_iterator end() const;
  std::size_t size() const;
  BlockMap_ptr_type ComponentMap(const std::string& name, bool ghosted=false) const;

  Comm_ptr_type Comm() const;

  bool HasComponent(const std::string& name) const;
  std::size_t NumComponents() const { return size(); }
  std::size_t NumVectors(const std::string& name) const;

  LO MyLength(const std::string& name, bool ghosted=false) const;
  GO GlobalLength(bool ghosted=false) const;
  
  // -- View data. --
  // Access a view of a single component's data.
  //
  // Const access -- this does not tag as changed.
  template<class DeviceType=AmanziDefaultDevice>
  cMultiVectorView_type_<DeviceType,Scalar>
  ViewComponent(const std::string& name, bool ghosted=false) const {
    return ghosted ? ghostvec_->ViewComponent(name) : mastervec_->ViewComponent(name);
  }

  //
  // Access a slice of a single component's data, at a specific dof num.
  //
  // Const access.
  template<class DeviceType=AmanziDefaultDevice>
  cVectorView_type_<DeviceType,Scalar>
  ViewComponent(const std::string& name, int dof, bool ghosted=false) const {
    return ghosted ? ghostvec_->ViewComponent(name, dof) : mastervec_->ViewComponent(name, dof);
  }

  // -- Set data. --

  // Access a view of a single component's data.
  //
  // Non-const access -- tags changed.
  template<class DeviceType=AmanziDefaultDevice>
  MultiVectorView_type_<DeviceType,Scalar>
  ViewComponent(const std::string& name, bool ghosted=false) {
    return ghosted ? ghostvec_->ViewComponent(name) : mastervec_->ViewComponent(name);
  }

  //
  // Access a slice of a single component's data, at a specific dof num.
  //
  // Non-const access.
  template<class DeviceType=AmanziDefaultDevice>
  VectorView_type_<DeviceType,Scalar>
  ViewComponent(const std::string& name, int dof, bool ghosted=false) {
    return ghosted ? ghostvec_->ViewComponent(name, dof) : mastervec_->ViewComponent(name, dof);
  }

  // Scatter master values to ghosted values, on all components, in a mode.
  //
  // Insert overwrites the current ghost value with the (unique) new
  // master value.
  //
  // Note that although scatter changes things, it doesn't change master
  // data, so we allow it to work on const.  This is necessary for a
  // non-owning PK to communicate a non-owned vector.
  //
  // This Scatter() is not managed, and is always done.  Tags changed.
  void ScatterMasterToGhosted(Tpetra::CombineMode mode=Tpetra::INSERT) const;
  void ScatterMasterToGhosted(const std::string& name, Tpetra::CombineMode=Tpetra::INSERT) const;

  // Combine ghosted values back to master values.
  //
  // Modes shown in Tpetra::CombineMode.h, but the default is ADD,
  // where off-process values are first summed into the on-process value.
  //
  // This Scatter() is not managed, and is always done.  Tags changed.
  void GatherGhostedToMaster(Tpetra::CombineMode mode=Tpetra::ADD);
  void GatherGhostedToMaster(const std::string& name, Tpetra::CombineMode mode=Tpetra::ADD);

  // returns non-empty importer
  Import_ptr_type importer(const std::string& name);

  // -- Assorted vector operations, this implements a Vec --

  // Sets all vectors to value.
  int PutScalar(Scalar scalar);

  // Sets all vectors to value including ghosted elements.
  // Different name is given so it cannot be used in a templated code.   
  int PutScalarMasterAndGhosted(Scalar scalar);

  // Sets ghost elements to value.
  // Different name is given so it cannot be used in a templated code.   
  int PutScalarGhosted(Scalar scalar);

  // v(name,:,:) = scalar
  int PutScalar(const std::string& name, Scalar scalar);

  // v(name,i,:) = scalar[i]
  int PutScalar(const std::string& name, std::vector<Scalar> scalar);

  // this <- scalar*this
  int Scale(Scalar scalar);
  int ScaleMasterAndGhosted(Scalar scalar);

  // this <- abs(this)
  int Abs(const CompositeVector_<Scalar>& other);
  
  // this(name,:,:) <- scalar*this(name,:,:)
  int Scale(const std::string& name, Scalar scalar);

  // // this <- this + scalar
  // int Shift(Scalar scalar);

  // // this(name,:,:) <- scalar + this(name,:,:)
  // int Shift(const std::string& name, Scalar scalar);

  // this <- element wise reciprocal(this)
  int Reciprocal(const CompositeVector_<Scalar>& other);
  
  // result <- other \dot this
  int Dot(const CompositeVector_<Scalar>& other, Scalar* result) const;

  // this <- scalarA*A + scalarThis*this
  CompositeVector_<Scalar>& Update(Scalar scalarA, const CompositeVector_<Scalar>& A, Scalar scalarThis);

  // this <- scalarA*A + scalarB*B + scalarThis*this
  CompositeVector_<Scalar>& Update(Scalar scalarA, const CompositeVector_<Scalar>& A,
                          Scalar scalarB, const CompositeVector_<Scalar>& B, Scalar scalarThis);

  // this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
  int Multiply(Scalar scalarAB, const CompositeVector_<Scalar>& A, const CompositeVector_<Scalar>& B,
               Scalar scalarThis);

  // this <- scalarAB * A^-1@B + scalarThis*this  (@ is the elementwise product
  // int ReciprocalMultiply(Scalar scalarAB, const CompositeVector_<Scalar>& A,
  //                        const CompositeVector_<Scalar>& B, Scalar scalarThis);

  // -- norms --
  int NormInf(Scalar* norm) const;
  int Norm1(Scalar* norm) const;
  int Norm2(Scalar* norm) const;

  // int MinValue(Scalar* value) const;
  // int MaxValue(Scalar* value) const;
  // int MeanValue(Scalar* value) const;

  // void MinValue(std::map<std::string, Scalar>& value) const;
  // void MaxValue(std::map<std::string, Scalar>& value) const;
  // void MeanValue(std::map<std::string, Scalar>& value) const;

  // -- Utilities --

  // Write components to outstream.
  void Print(std::ostream &os, bool data_io = true) const;

  // Populate by random numbers between -1 and 1.
  int Random();

 private:
  void InitMap_(const CompositeMap& space);
  void InitData_(const CompositeVector_<Scalar>& other, InitMode mode);
  void CreateData_();

  int Index_(const std::string& name) const {
    std::map<std::string, int>::const_iterator item = indexmap_.find(name);
    AMANZI_ASSERT(item != indexmap_.end());
    return item->second;
  }


  // The Vandelay is an Importer/Exporter which allows face unknowns
  // to be spoofed as boundary face unknowns.
  void CreateVandelay_() const;
  void ApplyVandelay_() const;

  // Tpetra vector accessors
  const MultiVector_type_<Scalar>&
  GetComponent_(const std::string& name, bool ghosted=false) const;
  MultiVector_type_<Scalar>&
  GetComponent_(const std::string& name, bool ghosted=false);

 private:
  Teuchos::RCP<const CompositeMap> map_;
  bool ghosted_;

  // data enumerating the blocks
  std::map< std::string, int > indexmap_;
  std::vector<std::string> names_;

  // data containing the blocks
  mutable Teuchos::RCP<CompVector<Scalar> > ghostvec_;
  Teuchos::RCP<CompVector<Scalar> > mastervec_;

  // importers for scatter/gather operation
  mutable std::vector<Import_ptr_type> importers_;

  // importer and vector for boundary data
  mutable Import_ptr_type vandelay_importer_;
  mutable MultiVector_ptr_type_<Scalar> vandelay_vector_;
};

} // namespace


#endif
