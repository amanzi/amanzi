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

} // namespace
