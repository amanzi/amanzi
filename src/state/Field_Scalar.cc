/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for a Field.  Field is not intended so much to hide implementation
of data as to restrict write access to data.  It freely passes out pointers to
its private data, but only passes out read-only const pointers unless you have
the secret password (AKA the name of the process kernel that owns the data).

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#include "errors.hh"
#include "Field.hh"
#include "Field_Scalar.hh"

namespace Amanzi {

Field_Scalar::Field_Scalar(std::string fieldname, std::string owner) :
    Field::Field(fieldname, owner), data_() {
  type_ = CONSTANT_SCALAR;
};

Field_Scalar::Field_Scalar(std::string fieldname, std::string owner,
                           Teuchos::RCP<double>& data) :
    Field::Field(fieldname, owner), data_(data) {
  type_ = CONSTANT_SCALAR;
};

// copy constructor:
Field_Scalar::Field_Scalar(const Field_Scalar& other) :
    Field::Field(other.fieldname_, other.owner_) {
  io_restart_ = other.io_restart_;
  io_vis_ = other.io_vis_;
  initialized_ = other.initialized_;
  type_ = other.type_;
  data_ = Teuchos::rcp(new double(*other.data_));
};

// Virtual copy constructor
Teuchos::RCP<Field> Field_Scalar::Clone() const {
  return Teuchos::rcp(new Field_Scalar(*this));
}

// Virtual copy constructor with non-empty name
Teuchos::RCP<Field> Field_Scalar::Clone(std::string fieldname) const {
  Teuchos::RCP<Field_Scalar> other = Teuchos::rcp(new Field_Scalar(*this));
  other->fieldname_ = fieldname;
  return other;
};

// Virtual copy constructor with non-empty name and new owner
Teuchos::RCP<Field> Field_Scalar::Clone(std::string fieldname, std::string owner)  const {
  Teuchos::RCP<Field_Scalar> other = Teuchos::rcp(new Field_Scalar(*this));
  other->fieldname_ = fieldname;
  other->owner_ = owner;
  return other;
};

// write-access to the data
Teuchos::RCP<double> Field_Scalar::scalar_data(std::string pk_name) {
  assert_owner_or_die_(pk_name);
  return data_;
};

// Overwrite data by pointer, not copy
void Field_Scalar::set_data(std::string pk_name, Teuchos::RCP<double>& data) {
  assert_owner_or_die_(pk_name);
  data_ = data;
};

// Set data by copy.
void Field_Scalar::set_data(std::string pk_name, double data) {
  assert_owner_or_die_(pk_name);
  *data_ = data;
};

// Initialization
void Field_Scalar::Initialize(Teuchos::ParameterList& plist) {
  if (plist.isParameter("Constant "+fieldname_)) {
    set_data(owner_, plist.get<double>("Constant "+(fieldname_)));
    set_initialized();
  }
};

} // namespace
