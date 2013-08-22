/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for a Field.  A field contains meta-data about a data structure and
that data structure, encapsulating things like vis, checkpointing, ownership,
initialization, etc.

------------------------------------------------------------------------- */

#ifndef AMANZI_STATE_FIELD_HH_
#define AMANZI_STATE_FIELD_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

#include "visualization.hh"
#include "checkpoint.hh"
#include "composite_vector.hh"

#include "state_defs.hh"

namespace Amanzi {

class Field {

public:
  // virtual destructor
  virtual ~Field() {}

  // copy constructor and assignment
  virtual Teuchos::RCP<Field> Clone() const = 0;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname) const = 0;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname, std::string owner) const = 0;

  // access
  std::string fieldname() const { return fieldname_; }
  std::string owner() const { return owner_; }
  std::string vis_key() const { return vis_key_; }
  FieldType type() const { return type_; }
  bool io_checkpoint() const { return io_checkpoint_; }
  bool io_vis() const { return io_vis_; }
  bool initialized() const { return initialized_; }

  // mutators
  void set_owner(std::string owner) { owner_ = owner; }
  void set_vis_key(std::string key) { vis_key_ = key; }
  void set_io_checkpoint(bool io_checkpoint=true) { io_checkpoint_ = io_checkpoint; }
  void set_io_vis(bool io_vis=true) { io_vis_ = io_vis; }
  void set_initialized(bool initialized=true) { initialized_ = initialized; }

  // Creation of the actual data (data is created lazily, allowing empty fields).
  virtual void CreateData() = 0;

  // data access and mutators
  virtual Teuchos::RCP<const double> GetScalarData() const {
    assert_type_or_die_(CONSTANT_SCALAR);
    not_implemented_error_();
  }
  virtual Teuchos::RCP<double> GetScalarData() {
    assert_type_or_die_(CONSTANT_SCALAR);
    not_implemented_error_();
  }
  virtual Teuchos::RCP<const Epetra_Vector> GetConstantVectorData() const {
    assert_type_or_die_(CONSTANT_VECTOR);
    not_implemented_error_();
  }
  virtual Teuchos::RCP<Epetra_Vector> GetConstantVectorData() {
    assert_type_or_die_(CONSTANT_VECTOR);
    not_implemented_error_();
  }
  virtual Teuchos::RCP<const CompositeVector> GetFieldData() const {
    assert_type_or_die_(COMPOSITE_VECTOR_FIELD);
    not_implemented_error_();
  }
  virtual Teuchos::RCP<CompositeVector> GetFieldData() {
    assert_type_or_die_(COMPOSITE_VECTOR_FIELD);
    not_implemented_error_();
  }

  // set data by pointer -- does not copy
  virtual void SetData(const Teuchos::RCP<CompositeVector>& data) {
    assert_type_or_die_(COMPOSITE_VECTOR_FIELD);
    not_implemented_error_();
  }
  virtual void SetData(const Teuchos::RCP<Epetra_Vector>& data) {
    assert_type_or_die_(CONSTANT_VECTOR);
    not_implemented_error_();
  }
  virtual void SetData(const Teuchos::RCP<double>& data) {
    assert_type_or_die_(CONSTANT_SCALAR);
    not_implemented_error_();
  }

  // set data by reference -- copies
  virtual void SetData(const CompositeVector& data) {
    assert_type_or_die_(COMPOSITE_VECTOR_FIELD);
    not_implemented_error_();
  }
  virtual void SetData(const Epetra_Vector& data) {
    assert_type_or_die_(CONSTANT_VECTOR);
    not_implemented_error_();
  }
  virtual void SetData(const double& data) {
    assert_type_or_die_(CONSTANT_SCALAR);
    not_implemented_error_();
  }

  // Initialize from a parameter list.
  virtual void Initialize(Teuchos::ParameterList& plist) = 0;

  // Initialize from a checkpoint file
  virtual void ReadCheckpoint(const Teuchos::Ptr<HDF5_MPI>& file) {}

  // write data to the visualization file
  virtual void WriteVis(const Teuchos::Ptr<Visualization>& vis) = 0;
  virtual void WriteCheckpoint(const Teuchos::Ptr<Checkpoint>& chk) = 0;

  // Compute from a function, if implemented by a subclass
  virtual void Compute(double time) {  }

protected:
  // constructors protected, should only be called by derived classes
  Field(std::string fieldname, std::string owner="state");
  Field(const Field& other);


  // consistency checking
  void assert_type_or_die_(FieldType type) const;
  void not_implemented_error_() const;

  FieldType type_;
  std::string fieldname_;
  std::string owner_;
  std::string vis_key_;

  bool io_checkpoint_;
  bool io_vis_;
  bool initialized_;

private:
  // operator= disabled
  Field& operator=(const Field&);

}; // class Field

} // namespace Amanzi

#endif
