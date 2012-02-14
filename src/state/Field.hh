/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for a Field.  Field is not intended so much to hide implementation
of data as to restrict write access to data.  It freely passes out pointers to
its private data, but only passes out read-only const pointers unless you have
the secret password (AKA the name of the process kernel that owns the data).

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#ifndef STATE_FIELD_HH_
#define STATE_FIELD_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

#include "Vis.hpp"
#include "CompositeVector.hh"

namespace Amanzi {

typedef enum { VECTOR_FIELD, CONSTANT_VECTOR, CONSTANT_SCALAR } FieldType;

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
  FieldType type() const { return type_; }
  bool io_restart() const { return io_restart_; }
  bool io_vis() const { return io_vis_; }
  bool initialized() const { return initialized_; }

  // mutators
  void set_owner(std::string owner) { owner_ = owner; }
  void set_io_restart(bool io_restart=true) { io_restart_ = io_restart; }
  void set_io_vis(bool io_vis=true) { io_vis_ = io_vis; }
  void set_initialized(bool initialized=true) { initialized_ = initialized; }

  // data access and mutators
  virtual Teuchos::RCP<const double> GetScalarData() const {
    assert_type_or_die_(CONSTANT_SCALAR);
    not_implemented_error_();
  }
  virtual Teuchos::RCP<double> GetScalarData(std::string pk_name) {
    assert_type_or_die_(CONSTANT_SCALAR);
    not_implemented_error_();
  }
  virtual Teuchos::RCP<const Epetra_Vector> GetConstantVectorData() const {
    assert_type_or_die_(CONSTANT_VECTOR);
    not_implemented_error_();
  }
  virtual Teuchos::RCP<Epetra_Vector> GetConstantVectorData(std::string pk_name) {
    assert_type_or_die_(CONSTANT_VECTOR);
    not_implemented_error_();
  }
  virtual Teuchos::RCP<const CompositeVector> GetFieldData() const {
    assert_type_or_die_(VECTOR_FIELD);
    not_implemented_error_();
  }
  virtual Teuchos::RCP<CompositeVector> GetFieldData(std::string pk_name) {
    assert_type_or_die_(VECTOR_FIELD);
    not_implemented_error_();
  }

  // set data by pointer -- does not copy
  virtual void SetData(std::string pk_name, Teuchos::RCP<CompositeVector>& data) {
    assert_type_or_die_(VECTOR_FIELD);
    not_implemented_error_();
  }
  virtual void SetData(std::string pk_name, Teuchos::RCP<Epetra_Vector>& data) {
    assert_type_or_die_(CONSTANT_VECTOR);
    not_implemented_error_();
  }
  virtual void SetData(std::string pk_name, Teuchos::RCP<double>& data) {
    assert_type_or_die_(CONSTANT_SCALAR);
    not_implemented_error_();
  }

  // set data by reference -- copies
  virtual void SetData(std::string pk_name, const CompositeVector& data) {
    assert_type_or_die_(VECTOR_FIELD);
    not_implemented_error_();
  }
  virtual void SetData(std::string pk_name, const Epetra_Vector& data) {
    assert_type_or_die_(CONSTANT_VECTOR);
    not_implemented_error_();
  }
  virtual void SetData(std::string pk_name, const double& data) {
    assert_type_or_die_(CONSTANT_SCALAR);
    not_implemented_error_();
  }

  // Initialize from a parameter list.
  virtual void Initialize(Teuchos::ParameterList& plist) = 0;

  // vis
  virtual void WriteVis(Amanzi::Vis& vis) = 0;

protected:
  // constructor protected, should only be called by derived classes
  Field(std::string fieldname, std::string owner="state");


  // consistency checking
  void assert_owner_or_die_(std::string pk_name) const;
  void assert_type_or_die_(FieldType type) const;
  void not_implemented_error_() const;

  FieldType type_;
  std::string fieldname_;
  std::string owner_;

  bool io_restart_;
  bool io_vis_;
  bool initialized_;

private:
  // operator= disabled
  Field& operator=(const Field&);

}; // class Field

} // namespace Amanzi

#endif
