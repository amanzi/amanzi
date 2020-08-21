/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for a single vector Field that is constant on the mesh.

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#include "errors.hh"
#include "Epetra_SerialComm.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"

#include "Field.hh"
#include "Field_ConstantVector.hh"

namespace Amanzi {

Field_ConstantVector::Field_ConstantVector(std::string fieldname, std::string owner) :
    Field::Field(fieldname, owner), data_(), dimension_(-1) {
  type_ = CONSTANT_VECTOR;
};

Field_ConstantVector::Field_ConstantVector(std::string fieldname, std::string owner,
        int dimension) :
    Field::Field(fieldname, owner), data_(), dimension_(dimension) {
  type_ = CONSTANT_VECTOR;
};

Field_ConstantVector::Field_ConstantVector(std::string fieldname, std::string owner,
                           Teuchos::RCP<Epetra_Vector>& data) :
    Field::Field(fieldname, owner), data_(data) {
  type_ = CONSTANT_VECTOR;
  dimension_ = data_->MyLength();
};

// copy constructor:
Field_ConstantVector::Field_ConstantVector(const Field_ConstantVector& other) :
    Field::Field(other),
    dimension_(other.dimension_),
    subfield_names_(other.subfield_names_) {
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


// meta data accesor
void Field_ConstantVector::set_dimension(int dimension) {
  if (dimension_ > 0 && dimension_ != dimension) {
      std::stringstream messagestream;
      messagestream << "In field " << fieldname_ <<
          " there were requests for sizes " << dimension_ << " and " << dimension << ".";
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
  }
  dimension_ = dimension;
}

// create data
void Field_ConstantVector::CreateData() {
  AMANZI_ASSERT(dimension_ > 0);
  Epetra_SerialComm comm;
  Epetra_LocalMap map(dimension_, 0, comm);
  data_ = Teuchos::rcp(new Epetra_Vector(map));
};

// write-access to the data
Teuchos::RCP<Epetra_Vector> Field_ConstantVector::GetConstantVectorData() {
  return data_;
};

// Overwrite data by pointer, not copy
void Field_ConstantVector::SetData(const Teuchos::RCP<Epetra_Vector>& data) {
  data_ = data;
};

// Set data by copy.
void Field_ConstantVector::SetData(const Epetra_Vector& data) {
  *data_ = data;
};


// Initialization
void Field_ConstantVector::Initialize(Teuchos::ParameterList& plist) {
  Teuchos::Array<double> vals = plist.get<Teuchos::Array<double> >("value");
  AMANZI_ASSERT(vals.size() == data_->MyLength());
  for (int i=0; i!=vals.size(); ++i) {
    (*data_)[i] = vals[i];
  }
  set_initialized();
};

void Field_ConstantVector::WriteVis(Visualization& vis) {}
void Field_ConstantVector::WriteCheckpoint(Checkpoint& chk) {}

} // namespace
