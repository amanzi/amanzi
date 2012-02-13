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
#include "Epetra_Vector.h"
#include "Field.hh"
#include "Field_ConstantVector.hh"

namespace Amanzi {

Field_ConstantVector::Field_ConstantVector(std::string fieldname, std::string owner) :
    Field::Field(fieldname, owner), data_() {
  type_ = CONSTANT_VECTOR;
};

Field_ConstantVector::Field_ConstantVector(std::string fieldname, std::string owner,
                           Teuchos::RCP<Epetra_Vector>& data) :
    Field::Field(fieldname, owner), data_(data) {
  type_ = CONSTANT_VECTOR;
};

// copy constructor:
Field_ConstantVector::Field_ConstantVector(const Field_ConstantVector& other) :
    Field::Field(other.fieldname_, other.owner_) {
  io_restart_ = other.io_restart_;
  io_vis_ = other.io_vis_;
  initialized_ = other.initialized_;
  type_ = other.type_;
  data_ = Teuchos::rcp(new Epetra_Vector(*other.data_));
};

// Virtual copy constructor
Teuchos::RCP<Field> Field_ConstantVector::Clone() const {
  return Teuchos::rcp(new Field_ConstantVector(*this));
}

// Virtual copy constructor with non-empty name
Teuchos::RCP<Field> Field_ConstantVector::Clone(std::string fieldname) const {
  Teuchos::RCP<Field_ConstantVector> other = Teuchos::rcp(new Field_ConstantVector(*this));
  other->fieldname_ = fieldname;
  return other;
};

// Virtual copy constructor with non-empty name and new owner
Teuchos::RCP<Field> Field_ConstantVector::Clone(std::string fieldname,
        std::string owner)  const {
  Teuchos::RCP<Field_ConstantVector> other = Teuchos::rcp(new Field_ConstantVector(*this));
  other->fieldname_ = fieldname;
  other->owner_ = owner;
  return other;
};

// write-access to the data
Teuchos::RCP<Epetra_Vector> Field_ConstantVector::constant_vector_data(std::string pk_name) {
  assert_owner_or_die_(pk_name);
  return data_;
};

// Overwrite data by pointer, not copy
void Field_ConstantVector::set_data(std::string pk_name,
        Teuchos::RCP<Epetra_Vector>& data) {
  assert_owner_or_die_(pk_name);
  data_ = data;
};

// Set data by copy.
void Field_ConstantVector::set_data(std::string pk_name, const Epetra_Vector& data) {
  assert_owner_or_die_(pk_name);
  *data_ = data;
};


// Initialization
void Field_ConstantVector::Initialize(Teuchos::ParameterList& plist) {
  if (subfield_names_.size() == data_->MyLength()) {
    Epetra_Vector vals(*data_);
    bool got_them_all = true;
    for (unsigned int lcv=0; lcv != subfield_names_.size(); ++lcv) {
      // attempt to pick out constant values for the field
      if (plist.isParameter("Constant "+subfield_names_[lcv])) {
        vals[lcv] = plist.get<double>("Constant "+subfield_names_[lcv]);
        std::cout << "  got value:" << fieldname_ << "=" << vals[lcv] << std::endl;
      } else {
        got_them_all = false;
        break;
      }
    }

    if (got_them_all) {
      std::cout << "got them all, assigning" << std::endl;
      *data_ = vals;
      set_initialized();
    }
  }
};

} // namespace
