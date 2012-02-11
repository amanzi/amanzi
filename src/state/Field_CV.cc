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
#include "CompositeVector.hh"
#include "Field.hh"
#include "Field_CV.hh"

namespace Amanzi {

Field_CV::Field_CV(std::string fieldname, std::string owner) :
    Field::Field(fieldname, owner), data_() {
  type_ = VECTOR_FIELD;
};

Field_CV::Field_CV(std::string fieldname, std::string owner,
                           Teuchos::RCP<CompositeVector>& data) :
    Field::Field(fieldname, owner), data_(data) {
  type_ = VECTOR_FIELD;
};

// copy constructor:
Field_CV::Field_CV(const Field_CV& other) :
    Field::Field(other.fieldname_, other.owner_) {
  io_restart_ = other.io_restart_;
  io_vis_ = other.io_vis_;
  initialized_ = other.initialized_;
  type_ = other.type_;
  data_ = Teuchos::rcp(new CompositeVector(*other.data_));
};

// Virtual copy constructor
Teuchos::RCP<Field> Field_CV::Clone() const {
  return Teuchos::rcp(new Field_CV(*this));
}

// Virtual copy constructor with non-empty name
Teuchos::RCP<Field> Field_CV::Clone(std::string fieldname) const {
  Teuchos::RCP<Field_CV> other = Teuchos::rcp(new Field_CV(*this));
  other->fieldname_ = fieldname;
  return other;
};

// Virtual copy constructor with non-empty name
Teuchos::RCP<Field> Field_CV::Clone(std::string fieldname, std::string owner) const {
  Teuchos::RCP<Field_CV> other = Teuchos::rcp(new Field_CV(*this));
  other->fieldname_ = fieldname;
  other->owner_ = owner;
  return other;
};

// write-access to the data
Teuchos::RCP<CompositeVector> Field_CV::vector_data(std::string pk_name) {
  assert_owner_or_die_(pk_name);
  return data_;
};

// Overwrite data by pointer, not copy
void Field_CV::set_data(std::string pk_name, Teuchos::RCP<CompositeVector>& data) {
  assert_owner_or_die_(pk_name);
  data_ = data;
};

void Field_CV::set_data(std::string pk_name, const CompositeVector& data) {
  assert_owner_or_die_(pk_name);
  *data_ = data;
};

} // namespace
