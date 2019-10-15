/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include "Record.hh"
#include "errors.hh"

namespace Amanzi {

// Basic constructor
Record::Record(Key fieldname, Key owner)
  : fieldname_(std::move(fieldname)),
    owner_(std::move(owner)),
    vis_key_(fieldname_),
    units_(),
    io_checkpoint_(true),
    io_vis_(true),
    initialized_(false)
{}

// Copy constructor does not copy data!
Record::Record(const Record& other)
  : fieldname_(other.fieldname_),
    owner_(other.owner_),
    vis_key_(other.vis_key_),
    units_(other.units_),
    io_checkpoint_(other.io_checkpoint_),
    io_vis_(other.io_vis_),
    initialized_(other.initialized_),
    subfieldnames_(other.subfieldnames_)
{}

// Record&
// Record::operator=(const Record& other)
// {
//   if (&other != this) {
//     fieldname_ = other.fieldname_;
//     owner_ = other.owner_;
//     vis_key_ = other.vis_key_;
//     units_ = other.units_;
//     io_checkpoint_ = other.io_checkpoint_;
//     io_vis_ = other.io_vis_;
//     initialized_ = other.initialized_;
//     subfieldnames_ = other.subfieldnames_;
//     data_ = other.data_;
//   }
// }

// pass-throughs for other functionality
void
Record::WriteVis(const Visualization& vis, Teuchos::ParameterList attrs) const
{
  if (io_vis()) {
    attributes(attrs);
    data_.WriteVis(vis, attrs);
  }
}
void
Record::WriteCheckpoint(const Checkpoint& chkp,
                        Teuchos::ParameterList attrs) const
{
  if (io_checkpoint()) {
    attributes(attrs);
    data_.WriteCheckpoint(chkp, attrs);
  }
}
void
Record::ReadCheckpoint(const Checkpoint& chkp, Teuchos::ParameterList attrs)
{
  if (io_checkpoint()) {
    attributes(attrs);
    data_.ReadCheckpoint(chkp, attrs);
  }
}

bool
Record::Initialize(Teuchos::ParameterList& plist, Teuchos::ParameterList attrs)
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
    attributes(attrs);
    initialized_here = data_.Initialize(plist, attrs);
    if (initialized_here) set_initialized();
  }
  return initialized_here;
}


void
Record::AssertOwnerOrDie(const Key& test_owner) const
{
  if (test_owner != owner()) {
    Errors::Message message;
    message << "Record \"" << fieldname_ << "\" requested by \"" << test_owner
            << "\" but is owned by \"" << owner() << "\"";
    throw(message);
  }
};

} // namespace Amanzi
