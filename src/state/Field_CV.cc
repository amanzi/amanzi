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

namespace Amanzi {

Field::Field(std::string fieldname, std::string owner) :
    fieldname_(fieldname),
    owner_(owner),
    io_restart_(true),
    io_vis_(false),
    initialized_(false),
    data_(NULL) {};

Field::Field(std::string fieldname, std::string owner, Teuchos::RCP<CompositeVector>& data):
    fieldname_(fieldname),
    owner_(owner),
    io_restart_(true),
    io_vis_(false),
    initialized_(false),
    data_(data) {};

// copy constructor:
// Create a new Field with different data but the same values.
Field::Field(const Field& field) {
  data_ = Teuchos::rcp(new CompositeVector(*field.data_));
  fieldname_ = field.fieldname_;
  owner_ = field.owner_;
  subfieldnames_ = field.subfieldnames_;
  io_restart_ = field.io_restart_;
  io_vis_ = field.io_vis_;
  initialized_ = field.initialized_;
};

// assignment is crucial... copy values not data/pointers
Field& Field::operator=(const Field& field) {
  if (this != &field) {
    *data_ = *(field.data_);
    fieldname_ = field.fieldname_;
    owner_ = field.owner_;
    subfieldnames_ = field.subfieldnames_;
    io_restart_ = field.io_restart_;
    io_vis_ = field.io_vis_;
    initialized_ = field.initialized_;
  }
  return *this;
};

// I miss decorators, do they exist in C++?
// check that the requesting pk owns the data
void Field::assert_owner_or_die_(std::string pk_name) const {
  if (pk_name != owner_) {
    std::stringstream messagestream;
    messagestream << "PK " << pk_name << " is attempting to write to " << fieldname_
                  << " which is owned by " << owner_;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
};

// write-access to the data
Teuchos::RCP<CompositeVector> Field::data(std::string pk_name) {
  assert_owner_or_die_(pk_name);
  return data_;
};

// Overwrite data by pointer, not copy
void Field::set_data(std::string pk_name, Teuchos::RCP<CompositeVector>& data) {
  assert_owner_or_die_(pk_name);
  data_ = data;
};

// Overwrite data by copy.
void Field::set_data(std::string pk_name, CompositeVector& data) {
  assert_owner_or_die_(pk_name);
  *data_ = data;
};

} // namespace
