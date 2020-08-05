/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for a scalar Field.

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#include "errors.hh"
#include "Function.hh"
#include "FunctionFactory.hh"

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
    Field::Field(other),
    func_(other.func_) {
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


// data creation
void Field_Scalar::CreateData() {
  data_ = Teuchos::rcp(new double);
}


// write-access to the data
Teuchos::RCP<double> Field_Scalar::GetScalarData() {
  return data_;
};

// Overwrite data by pointer, not copy
void Field_Scalar::SetData(const Teuchos::RCP<double>& data) {
  data_ = data;
};

// Set data by copy.
void Field_Scalar::SetData(const double& data) {
  *data_ = data;
};

// Initialization
void Field_Scalar::Initialize(Teuchos::ParameterList& plist) {
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
void Field_Scalar::WriteVis(Visualization& vis) {};

// checkpoint
void Field_Scalar::WriteCheckpoint(Checkpoint& chk) {}

// Compute from a function
void Field_Scalar::Compute(double time) {
  if (func_ != Teuchos::null) {
    std::vector<double> x(1,time);
    SetData((*func_)(x));
  }
}

} // namespace
