/*
  Data Structures

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Interface for CompVector, an implementation of a slightly improved
  Epetra_MultiVector which spans multiple simplices and knows how to
  communicate itself.
*/

#ifndef AMANZI_COMPVECTOR_DECL_HH_
#define AMANZI_COMPVECTOR_DECL_HH_

#include <vector>
#include "Teuchos_RCP.hpp"
#include "AmanziTypes.hh"
#include "CompMap.hh"

namespace Amanzi {

template<typename Scalar>
class CompVector {

public:
  //
  // Constructors etc
  // --------------------------------
  CompVector(const Teuchos::RCP<const CompMap>& map);

  // copy constructor
  CompVector(const CompVector& other);

  // Constructor just does maps, this creates data.
  void CreateData();

  // assignment
  CompVector& operator=(const CompVector<Scalar>& other);

  //
  // Accessors
  // --------------------------------
  Comm_ptr_type Comm() const { return Map()->Comm(); }

  // Iteration over names of the vector
  using name_iterator = CompMap::name_iterator;
  name_iterator begin() { return Map()->begin(); }
  name_iterator end() { return Map()->end(); }
  std::size_t size() const { return Map()->size(); }
  LO MyLength(const std::string& name) const;

  GO GlobalLength() const { return Map()->GlobalLength(); }
  std::size_t NumComponents() const { return Map()->NumComponents(); }
  std::size_t NumVectors(const std::string& name) const { return Map()->NumVectors(name); }

  Teuchos::RCP<const CompMap> Map() const { return map_; }
  BlockMap_ptr_type ComponentMap(const std::string& name) const { return map_->ComponentMap(name); }

  // Accessors to data.
  bool HasComponent(const std::string& name) const { return map_->HasComponent(name); }

  // -- Access a component vector
  cMultiVector_ptr_type_<Scalar> GetComponent(const std::string& name) const;
  MultiVector_ptr_type_<Scalar> GetComponent(const std::string& name);

  // -- View a component vector
  template<class DeviceType=AmanziDefaultDevice>
  cMultiVectorView_type_<DeviceType,Scalar> ViewComponent(const std::string& name) const;
  template<class DeviceType=AmanziDefaultDevice>
  MultiVectorView_type_<DeviceType,Scalar> ViewComponent(const std::string& name);

  // -- SubView of a component vector
  template<class DeviceType=AmanziDefaultDevice>
  cVectorView_type_<DeviceType,Scalar> ViewComponent(const std::string& name, std::size_t dof) const;
  template<class DeviceType=AmanziDefaultDevice>
  VectorView_type_<DeviceType,Scalar> ViewComponent(const std::string& name, std::size_t dof);
  
  // Mutators of data
  // -- Set entries in the vectors.
  void SetComponent(const std::string& name, const MultiVector_ptr_type_<Scalar>& data);

  //
  // Vector operations.
  // --------------------------------
  // Insert value into data.
  int PutScalar(Scalar scalar);

  // Insert values into data, by DOF, not by component!
  //int PutScalar(std::vector<Scalar> scalar);

  // Insert value into component [name].
  int PutScalar(const std::string& name, Scalar scalar);

  // Insert values into component [name].
  int PutScalar(const std::string& name, const std::vector<Scalar>& scalar);

  // this <- abs(this)
  int Abs(const CompVector<Scalar>& other);

  // this <- this * scalarThis
  int Scale(Scalar value);

  // Scale() applied to component name.
  int Scale(const std::string& name, Scalar scalarThis);

  // // this <- this + scalarA
  // int Shift(Scalar scalarA);

  // // Shift() applied to component name.
  // int Shift(const std::string& name, Scalar scalarA);

  // this <- element wise reciprocal(this)
  int Reciprocal(const CompVector<Scalar>& other);
  
  // result <- other \dot this
  int Dot(const CompVector<Scalar>& other, Scalar* result) const;

  // this <- scalarA*A + scalarThis*this
  CompVector<Scalar>& Update(Scalar scalarA, const CompVector<Scalar>& A, Scalar scalarThis);

  // this <- scalarA*A + scalarB*B + scalarThis*this
  CompVector<Scalar>& Update(Scalar scalarA, const CompVector<Scalar>& A,
                          Scalar scalarB, const CompVector<Scalar>& B, Scalar scalarThis);

  // this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
  int Multiply(Scalar scalarAB, const CompVector<Scalar>& A, const CompVector<Scalar>& B,
               Scalar scalarThis);

  // this <- scalarAB * B/A + scalarThis*this  (/ is the elementwise division
  // int ReciprocalMultiply(Scalar scalarAB, const CompVector<Scalar>& A, const CompVector<Scalar>& B,
  //                        Scalar scalarThis);

  // Norms.
  int NormInf(Scalar* norm) const;
  int Norm1(Scalar* norm) const;
  int Norm2(Scalar* norm) const;

  // int MinValue(Scalar* value) const;
  // int MaxValue(Scalar* value) const;
  // int MeanValue(Scalar* value) const;

  // Extras
  void Print(std::ostream &os, bool data_io = true) const;

  int Random();

 protected:
  cMultiVector_ptr_type_<Scalar> GetComponent_(const std::string& name) const { return data_.at(name); }
  MultiVector_ptr_type_<Scalar> GetComponent_(const std::string& name) { return data_.at(name); }

 protected:
  Teuchos::RCP<const CompMap> map_;
  std::map<std::string, MultiVector_ptr_type_<Scalar> > data_;

};

} // namespace

#include "CompVector_impl.hh"

#endif
