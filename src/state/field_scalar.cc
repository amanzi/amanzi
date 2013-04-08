/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for a scalar Field.

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#include "errors.hh"
#include "function.hh"
#include "function-factory.hh"

#include "field.hh"
#include "field_scalar.hh"

namespace Amanzi {

FieldScalar::FieldScalar(std::string fieldname, std::string owner) :
    Field::Field(fieldname, owner), data_() {
  type_ = CONSTANT_SCALAR;
};

FieldScalar::FieldScalar(std::string fieldname, std::string owner,
                           Teuchos::RCP<double>& data) :
    Field::Field(fieldname, owner), data_(data) {
  type_ = CONSTANT_SCALAR;
};

// copy constructor:
FieldScalar::FieldScalar(const FieldScalar& other) :
    Field::Field(other),
    func_(other.func_) {
  data_ = Teuchos::rcp(new double(*other.data_));
};

// Virtual copy constructor
Teuchos::RCP<Field> FieldScalar::Clone() const {
  return Teuchos::rcp(new FieldScalar(*this));
}

// Virtual copy constructor with non-empty name
Teuchos::RCP<Field> FieldScalar::Clone(std::string fieldname) const {
  Teuchos::RCP<FieldScalar> other = Teuchos::rcp(new FieldScalar(*this));
  other->fieldname_ = fieldname;
  return other;
};

// Virtual copy constructor with non-empty name and new owner
Teuchos::RCP<Field> FieldScalar::Clone(std::string fieldname, std::string owner)  const {
  Teuchos::RCP<FieldScalar> other = Teuchos::rcp(new FieldScalar(*this));
  other->fieldname_ = fieldname;
  other->owner_ = owner;
  return other;
};


// data creation
void FieldScalar::CreateData() {
  data_ = Teuchos::rcp(new double);
}


// write-access to the data
Teuchos::RCP<double> FieldScalar::GetScalarData() {
  return data_;
};

// Overwrite data by pointer, not copy
void FieldScalar::SetData(const Teuchos::RCP<double>& data) {
  data_ = data;
};

// Set data by copy.
void FieldScalar::SetData(const double& data) {
  *data_ = data;
};

// Initialization
void FieldScalar::Initialize(Teuchos::ParameterList& plist) {
  if (plist.isParameter("value")) {
    SetData(plist.get<double>("value"));
    set_initialized();
  } else if (plist.isParameter("function")) {
    FunctionFactory fac;
    func_ = Teuchos::rcp(fac.Create(plist.sublist("function")));
    set_initialized();
  }
};

// visualization
void FieldScalar::WriteVis(const Teuchos::Ptr<Visualization>& vis) {
  if (io_vis_) {
  }
};

// checkpoint
void FieldScalar::WriteCheckpoint(const Teuchos::Ptr<Checkpoint>& chk) {
  if (io_checkpoint_) {
  }
};

// Compute from a function
void FieldScalar::Compute(double time) {
  if (func_ != Teuchos::null) {
    SetData((*func_)(&time));
  }
}

} // namespace
