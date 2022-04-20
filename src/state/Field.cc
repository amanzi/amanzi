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

// constructor
Field::Field(std::string fieldname, std::string owner) :
    type_(NULL_FIELD_TYPE),
    fieldname_(fieldname),
    owner_(owner),
    io_checkpoint_(true),
    io_vis_(true),
    initialized_(false),
    vis_key_("domain")
{};

Field::Field(const Field& other) :
    type_(other.type_),
    fieldname_(other.fieldname_),
    owner_(other.owner_),
    io_checkpoint_(other.io_checkpoint_),
    io_vis_(other.io_vis_),
    initialized_(other.initialized_),
    vis_key_(other.vis_key_)
{};

void Field::assert_type_or_die_(FieldType type) const {
  if (type != type_) {
    std::stringstream messagestream;
    messagestream << "Attempting to use a method intended for a type " << type
                  << " field on a field of type " << type_;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
};

void Field::not_implemented_error_() const {
  std::stringstream messagestream;
  messagestream << "Method not implemented for this type... blaim ecoon@lanl.gov!";
  Errors::Message message(messagestream.str());
  Exceptions::amanzi_throw(message);
};

bool Field::HasCopy(Key timetag) const {

  return GetCopy_(timetag) != Teuchos::null;

}

Teuchos::RCP<const Field> Field::GetCopy_(Key timetag) const{

  FieldMap::const_iterator lb = field_copy_.lower_bound(timetag);
  if (lb != field_copy_.end() && !(field_copy_.key_comp()(timetag, lb->first))) {
    return lb->second;
  } else {
    return Teuchos::null;
  }  
  
}

Teuchos::RCP<Field> Field::GetCopy_(Key timetag) {

  FieldMap::iterator lb = field_copy_.lower_bound(timetag);
  if (lb != field_copy_.end() && !(field_copy_.key_comp()(timetag, lb->first))) {
    return lb->second;
  } else {
  return Teuchos::null;
  }  
  
}


Teuchos::RCP<Field> Field::GetCopy(Key timetag, Key pk_name) {

  Teuchos::RCP<Field> record = GetCopy_(timetag);

  if (record == Teuchos::null) {
    std::stringstream messagestream;
    messagestream << "Field " << fieldname_ << " does not have copy tagged " << timetag;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  } else if (record->owner() != pk_name) {
    std::stringstream messagestream;
    messagestream << "PK \"" << pk_name
                  << "\" is attempting write access to copy \"" <<timetag
                  << "\" which is owned by \"" << GetCopy_(timetag)->owner() <<"\"";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  return record;

}

Teuchos::RCP<const Field>  Field::GetCopy(Key timetag) const {

  Teuchos::RCP<const Field> record = GetCopy_(timetag);
  
  if (record == Teuchos::null) {
    std::stringstream messagestream;
    messagestream << "Field " << fieldname_ << " does not have copy tagged " << timetag;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  return record;

}

// //  No data movement on pointer switch
// void Field::SwitchCopies(Key timetag1, Key timetag2) {

//   if (timetag1 != "default") {
//     if (!HasCopy(timetag1)) {
//       std::stringstream messagestream;
//       messagestream << "Field " << fieldname_ << " does not have copy tagged " << timetag1;
//       Errors::Message message(messagestream.str());
//       Exceptions::amanzi_throw(message);
//     }
//   }

//   if (timetag2 != "default") {
//     if (!HasCopy(timetag2)) {
//       std::stringstream messagestream;
//       messagestream << "Field " << fieldname_ << " does not have copy tagged " << timetag2;
//       Errors::Message message(messagestream.str());
//       Exceptions::amanzi_throw(message);
//     }
//   }


//   FieldMap::iterator lb1 = field_copy_.lower_bound(timetag1);
//   FieldMap::iterator lb2 = field_copy_.lower_bound(timetag2);

//   Teuchos::RCP<Field> record;
//   if (timetag1 != "default") {
//     record = lb1 ->second;
//   }else{
    
    
//   lb1->second = lb2->second;
//   lb2->second = record;

// }

// set data by pointer -- does not copy
void Field::SetCopy(Key timetag, const Teuchos::RCP<Field>& field) {

  assert_type_or_die_(field->type());

  if (!HasCopy(timetag)) {
    std::stringstream messagestream;
    messagestream << "Field " << fieldname_ << " does not have copy tagged " << timetag;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  FieldMap::iterator lb = field_copy_.lower_bound(timetag);

  lb->second = field;

}

void Field::RequireCopy(Key tag) {

  if (!HasCopy(tag)) {
    Teuchos::RCP<Field> copy = Clone(fieldname_, owner_);
    field_copy_[tag] = copy;
  }
}

void Field::RequireCopy(Key tag, Key new_owner) {

  if (!HasCopy(tag)) {
    Teuchos::RCP<Field> copy = Clone(fieldname_, new_owner);
    field_copy_[tag] = copy;
  }
}


} // namespace
