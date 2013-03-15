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

#include "field.hh"
#include "field_constant_vector.hh"

namespace Amanzi {

FieldConstantVector::FieldConstantVector(std::string fieldname, std::string owner) :
    Field::Field(fieldname, owner), data_(), dimension_(-1) {
  type_ = CONSTANT_VECTOR;
};

FieldConstantVector::FieldConstantVector(std::string fieldname, std::string owner,
        int dimension) :
    Field::Field(fieldname, owner), data_(), dimension_(dimension) {
  type_ = CONSTANT_VECTOR;
};

FieldConstantVector::FieldConstantVector(std::string fieldname, std::string owner,
                           Teuchos::RCP<Epetra_Vector>& data) :
    Field::Field(fieldname, owner), data_(data) {
  type_ = CONSTANT_VECTOR;
  dimension_ = data_->MyLength();
};

// copy constructor:
FieldConstantVector::FieldConstantVector(const FieldConstantVector& other) :
    Field::Field(other),
    dimension_(other.dimension_),
    subfield_names_(other.subfield_names_) {
  data_ = Teuchos::rcp(new Epetra_Vector(*other.data_));
};

// Virtual copy constructor
Teuchos::RCP<Field> FieldConstantVector::Clone() const {
  return Teuchos::rcp(new FieldConstantVector(*this));
}

// Virtual copy constructor with non-empty name
Teuchos::RCP<Field> FieldConstantVector::Clone(std::string fieldname) const {
  Teuchos::RCP<FieldConstantVector> other = Teuchos::rcp(new FieldConstantVector(*this));
  other->fieldname_ = fieldname;
  return other;
};

// Virtual copy constructor with non-empty name and new owner
Teuchos::RCP<Field> FieldConstantVector::Clone(std::string fieldname,
        std::string owner)  const {
  Teuchos::RCP<FieldConstantVector> other = Teuchos::rcp(new FieldConstantVector(*this));
  other->fieldname_ = fieldname;
  other->owner_ = owner;
  return other;
};


// meta data accesor
void FieldConstantVector::set_dimension(int dimension) {
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
void FieldConstantVector::CreateData() {
  ASSERT(dimension_ > 0);
  Epetra_SerialComm comm;
  Epetra_LocalMap map(dimension_, 0, comm);
  data_ = Teuchos::rcp(new Epetra_Vector(map));
};

// write-access to the data
Teuchos::RCP<Epetra_Vector> FieldConstantVector::GetConstantVectorData() {
  return data_;
};

// Overwrite data by pointer, not copy
void FieldConstantVector::SetData(const Teuchos::RCP<Epetra_Vector>& data) {
  data_ = data;
};

// Set data by copy.
void FieldConstantVector::SetData(const Epetra_Vector& data) {
  *data_ = data;
};


// Initialization
void FieldConstantVector::Initialize(Teuchos::ParameterList& plist) {
  Teuchos::Array<double> vals = plist.get<Teuchos::Array<double> >("value");
  ASSERT(vals.size() == data_->MyLength());
  for (int i=0; i!=vals.size(); ++i) {
    (*data_)[i] = vals[i];
  }
  set_initialized();
};

void FieldConstantVector::WriteVis(const Teuchos::Ptr<Visualization>& vis) {
  if (io_vis_) {
  }
};

void FieldConstantVector::WriteCheckpoint(const Teuchos::Ptr<Checkpoint>& chk) {
  if (io_checkpoint_) {
  }
};

} // namespace
