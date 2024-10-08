/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

  Interface for a Record. A field contains meta-data about a
  data structure such as ownership, initialization, etc.
*/

#include "Record.hh"
#include "errors.hh"

namespace Amanzi {

// Basic constructor
Record::Record(Key fieldname, Tag tag, Key owner)
  : fieldname_(std::move(fieldname)),
    tag_(std::move(tag)),
    owner_(std::move(owner)),
    vis_key_(fieldname_),
    initialized_(false),
    io_checkpoint_(true),
    io_vis_(true)
{}

// Copy constructor does not copy data!
Record::Record(const Record& other, const Tag* tag)
  : fieldname_(other.fieldname_),
    tag_(other.tag_),
    owner_(other.owner_),
    vis_key_(other.vis_key_),
    initialized_(other.initialized_),
    io_checkpoint_(other.io_checkpoint_),
    io_vis_(other.io_vis_)
{
  if (tag) tag_ = *tag;
}

// pass-throughs for other functionality
void
Record::WriteVis(const Visualization& vis, Teuchos::ParameterList& attrs) const
{
  if (io_vis()) { data_.WriteVis(vis, attrs); }
}

void
Record::WriteCheckpoint(const Checkpoint& chkp,
                        Teuchos::ParameterList& attrs,
                        bool post_mortem) const
{
  if (post_mortem || io_checkpoint()) data_.WriteCheckpoint(chkp, attrs);
}

bool
Record::ReadCheckpoint(const Checkpoint& chkp, Teuchos::ParameterList& attrs)
{
  if (io_checkpoint()) {
    data_.ReadCheckpoint(chkp, attrs);
    set_initialized();
    return true;
  }
  return false;
}

bool
Record::Initialize(Teuchos::ParameterList& plist, bool force)
{
  // check meta-data
  if (plist.isParameter("write checkpoint")) {
    set_io_checkpoint(plist.get<bool>("write checkpoint"));
  }

  if (plist.isParameter("write vis")) { set_io_vis(plist.get<bool>("write vis")); }

  bool initialized_here = false;
  if (!initialized() || force) {
    initialized_here = data_.Initialize(plist);
    if (initialized_here) set_initialized();
  }
  return initialized_here;
}


void
Record::Assign(const Record& other)
{
  try {
    data_.Assign(other.data_);
  } catch (const Errors::Message& msg) {
    Errors::Message new_msg;
    new_msg << "Access to field: \"" << fieldname() << "\"" << msg.what();
    throw(new_msg);
  }
}

void
Record::AssignPtr(const Record& other)
{
  try {
    data_ = other.data_;
  } catch (const Errors::Message& msg) {
    Errors::Message new_msg;
    new_msg << "Access to field: \"" << fieldname() << "\"" << msg.what();
    throw(new_msg);
  }
}

void
Record::AssertOwnerOrDie(const Key& test_owner) const
{
  if (test_owner != owner()) {
    Errors::Message message;
    message << "Record \"" << fieldname_ << "@" << tag_.get() << "\" requested by \"" << test_owner
            << "\" but is owned by \"" << owner() << "\"";
    throw(message);
  }
};

} // namespace Amanzi
