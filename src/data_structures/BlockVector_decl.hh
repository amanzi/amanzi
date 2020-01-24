/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

/*
  Interface for BlockVector, an implementation of a slightly improved
  Epetra_MultiVector which spans multiple simplices and knows how to
  communicate itself.
*/

#ifndef AMANZI_BLOCK_VECTOR_DECL_HH_
#define AMANZI_BLOCK_VECTOR_DECL_HH_

#include <vector>
#include "Teuchos_RCP.hpp"
#include "Teuchos_DataAccess.hpp"
#include "AmanziTypes.hh"
#include "DataStructuresHelpers.hh"
#include "BlockSpace.hh"

namespace Amanzi {

template <typename Scalar>
class BlockVector {
 public:
  //
  // Constructors etc
  // ---------------------------------------------
  BlockVector(const Teuchos::RCP<const BlockSpace>& map,
              InitMode mode = InitMode::ZERO);

  // copy constructor
  BlockVector(const BlockVector& other,
              Teuchos::DataAccess access = Teuchos::DataAccess::Copy,
              InitMode mode = InitMode::COPY);

  // assignment
  BlockVector<Scalar>& operator=(const BlockVector<Scalar>& other);
  void assign(const BlockVector<Scalar>& other) { *this = other; }

  //
  // Meta-data delegated to map
  // ---------------------------------------------
  const Teuchos::RCP<const BlockSpace>& getMap() const { return map_; }
  Comm_ptr_type Comm() const { return getMap()->Comm(); }
  GO getGlobalLength(bool ghosted = false) const
  {
    return getMap()->getGlobalLength(ghosted);
  }
  LO getLocalLength(bool ghosted = false) const
  {
    return getMap()->getLocalLength(ghosted);
  }

  //
  // Component meta-data delegated to map
  // ---------------------------------------------
  bool HasComponent(const std::string& name) const
  {
    return map_->HasComponent(name);
  }

  using name_iterator = BlockSpace::name_iterator;
  name_iterator begin() const { return getMap()->begin(); }
  name_iterator end() const { return getMap()->end(); }
  std::size_t size() const { return getMap()->size(); }

  std::size_t getNumVectors(const std::string& name) const
  {
    return getMap()->getNumVectors(name);
  }

  //
  // Accessors to data.
  // ---------------------------------------------
  template <class DeviceType = AmanziDefaultDevice>
  using cMultiVectorView_type = cMultiVectorView_type_<DeviceType, Scalar>;
  template <class DeviceType = AmanziDefaultDevice>
  using MultiVectorView_type = MultiVectorView_type_<DeviceType, Scalar>;
  template <class DeviceType = AmanziDefaultDevice>
  using cVectorView_type = cVectorView_type_<DeviceType, Scalar>;
  template <class DeviceType = AmanziDefaultDevice>
  using VectorView_type = VectorView_type_<DeviceType, Scalar>;

  // -- Access a component vector
  cMultiVector_ptr_type_<Scalar>
  GetComponent(const std::string& name, bool ghosted = false) const;
  MultiVector_ptr_type_<Scalar>
  GetComponent(const std::string& name, bool ghosted = false);

  // -- View a component vector
  template <class DeviceType = AmanziDefaultDevice>
  cMultiVectorView_type<DeviceType>
  ViewComponent(const std::string& name, bool ghosted = false) const;
  template <class DeviceType = AmanziDefaultDevice>
  MultiVectorView_type<DeviceType>
  ViewComponent(const std::string& name, bool ghosted = false);

  // -- SubView of a component vector
  template <class DeviceType = AmanziDefaultDevice>
  cVectorView_type<DeviceType>
  ViewComponent(const std::string& name, std::size_t dof,
                bool ghosted = false) const;
  template <class DeviceType = AmanziDefaultDevice>
  VectorView_type<DeviceType>
  ViewComponent(const std::string& name, std::size_t dof, bool ghosted = false);

  // // -- Set entries in the vectors.
  // void SetComponent(const std::string& name, const
  // MultiVector_ptr_type_<Scalar>& data);

  //
  // Communication operations
  // --------------------------------
  // Scatter master values to ghosted values, on all components, in a mode.
  //
  // Insert overwrites the current ghost value with the (unique) new
  // master value.
  //
  // Note that although scatter changes things, it doesn't change master
  // data, so we allow it to work on const.  This is necessary for a
  // non-owning PK to communicate a non-owned vector.
  void ScatterMasterToGhosted(Tpetra::CombineMode mode = Tpetra::INSERT) const;
  void ScatterMasterToGhosted(const std::string& name,
                              Tpetra::CombineMode = Tpetra::INSERT) const;

  // Combine ghosted values back to master values.
  //
  // Modes shown in Tpetra::CombineMode.h, but the default is ADD,
  // where off-process values are first summed into the on-process value.
  void GatherGhostedToMaster(Tpetra::CombineMode mode = Tpetra::ADD);
  void GatherGhostedToMaster(const std::string& name,
                             Tpetra::CombineMode mode = Tpetra::ADD);


  //
  // Vector operations.
  // --------------------------------
  // Insert value into data.
  void putScalar(Scalar scalar);

  // Insert value into component [name].
  void putScalar(const std::string& name, Scalar scalar);

  // Insert values into component [name].
  void putScalar(const std::string& name, const std::vector<Scalar>& scalar);

  // Sets all vectors to value including ghosted elements.
  void putScalarMasterAndGhosted(Scalar scalar);

  // Sets ghost elements to value.
  void putScalarGhosted(Scalar scalar);

  // cheap randomizer
  void random();

  // this <- abs(this)
  void abs(const BlockVector<Scalar>& other);

  // this <- this * scalarThis
  void scale(Scalar value);

  // scale() applied to component name.
  void scale(const std::string& name, Scalar scalarThis);

  // // this <- this + scalarA
  void shift(Scalar scalarA) {
    BlockVector<Scalar> one(*this);
    one.putScalar(1.0);
    this->update(scalarA, one, 1.);
  }

  // // Shift() applied to component name.
  // void Shift(const std::string& name, Scalar scalarA);

  // this <- element wise reciprocal(this)
  void reciprocal(const BlockVector<Scalar>& other);

  // result <- other \dot this
  Scalar dot(const BlockVector<Scalar>& other) const;

  // this <- scalarA*A + scalarThis*this
  void update(Scalar scalarA, const BlockVector<Scalar>& A, Scalar scalarThis);

  // this <- scalarA*A + scalarB*B + scalarThis*this
  void update(Scalar scalarA, const BlockVector<Scalar>& A, Scalar scalarB,
              const BlockVector<Scalar>& B, Scalar scalarThis);

  // this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
  void elementWiseMultiply(Scalar scalarAB, const BlockVector<Scalar>& A,
                           const BlockVector<Scalar>& B, Scalar scalarThis);

  // this <- scalarAB * B/A + scalarThis*this  (/ is the elementwise division
  // int ReciprocalelementWiseMultiply(Scalar scalarAB, const
  // BlockVector<Scalar>& A, const BlockVector<Scalar>& B,
  //                        Scalar scalarThis);

  // Norms.
  Scalar normInf() const;
  Scalar norm1() const;
  Scalar norm2() const;

  // Scalar MinValue() const;
  // Scalar MaxValue() const;
  // Scalar MeanValue() const;

  // Extras
  void Print(std::ostream& os, bool ghosted = false, bool data_io = true) const;

 protected:
  virtual cMultiVector_ptr_type_<Scalar>
  GetComponent_(const std::string& name, bool ghosted = false) const
  {
    return ghosted ? ghost_data_.at(name) : master_data_.at(name);
  }
  virtual MultiVector_ptr_type_<Scalar>
  GetComponent_(const std::string& name, bool ghosted = false)
  {
    return ghosted ? ghost_data_.at(name) : master_data_.at(name);
  }

  void SetComponent_(const std::string& name, bool ghosted,
                     const MultiVector_ptr_type_<Scalar>& v)
  {
    if (ghosted)
      ghost_data_[name] = v;
    else
      master_data_[name] = v;
  }

  // Constructor just does maps, this allocates memory.
  void CreateData_(InitMode mode);

 protected:
  Teuchos::RCP<const BlockSpace> map_;
  std::map<std::string, MultiVector_ptr_type_<Scalar>> master_data_,
    ghost_data_;
};

} // namespace Amanzi

#endif
