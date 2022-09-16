/*
  State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

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
     io_vis_(true) {}

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
void Record::WriteVis(const Visualization& vis,
                      const std::vector<std::string>* subfieldnames) const
{
  if (io_vis()) {
    data_.WriteVis(vis, vis_fieldname(), subfieldnames);
  }
}

void Record::WriteCheckpoint(const Checkpoint& chkp, const Tag& tag, bool post_mortem,
                             const std::vector<std::string>* subfieldnames) const
{
  if (post_mortem || io_checkpoint()) {
    auto name = Keys::getKey(vis_fieldname(), tag);
    data_.WriteCheckpoint(chkp, name, subfieldnames);
  }
}

bool Record::ReadCheckpoint(const Checkpoint& chkp, const Tag& tag,
                            const std::vector<std::string>* subfieldnames)
{
  if (io_checkpoint()) {
    auto name = Keys::getKey(vis_fieldname(), tag);
    data_.ReadCheckpoint(chkp, name, subfieldnames);
    return true;
  }
  return false;
}

bool Record::Initialize(Teuchos::ParameterList& plist,
                        const std::vector<std::string>* subfieldnames)
{
  // check meta-data
  if (plist.isParameter("write checkpoint")) {
    bool checkpoint_io = plist.get<bool>("write checkpoint", false);
    set_io_checkpoint(checkpoint_io);
  }

  if (plist.isParameter("write vis")) {
    bool vis_io = plist.get<bool>("write vis", false);
    set_io_vis(vis_io);
  }

  bool initialized_here = false;
  if (!initialized()) {
    initialized_here = data_.Initialize(plist, vis_fieldname(), subfieldnames);
    if (initialized_here) set_initialized();
  }
  return initialized_here;
}


void Record::Assign(const Record& other)
{
  try {
    data_.Assign(other.data_);
  } catch (const Errors::Message& msg) {
    Errors::Message new_msg;
    new_msg << "Access to field: \"" << fieldname() << "\"" << msg.what();
    throw(new_msg);
  }
}

void Record::AssignPtr(const Record& other)
{
  try {
    data_ = other.data_;
  } catch (const Errors::Message& msg) {
    Errors::Message new_msg;
    new_msg << "Access to field: \"" << fieldname() << "\"" << msg.what();
    throw(new_msg);
  }
}

void Record::AssertOwnerOrDie(const Key& test_owner) const
{
  if (test_owner != owner()) {
    Errors::Message message;
    message << "Record \"" << fieldname_ << "@" << tag_.get()
            << "\" requested by \"" << test_owner
            << "\" but is owned by \"" << owner() << "\"";
    throw(message);
  }
};

}  // namespace Amanzi
