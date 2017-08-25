/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon


  Interface for a Record.  A field contains meta-data about a data structure
  such as ownership, initialization, etc.
------------------------------------------------------------------------- */

#include "errors.hh"
#include "Record.hh"

namespace Amanzi {

// Basic constructor
Record::Record(Key fieldname, Key owner) :
    fieldname_(std::move(fieldname)),
    owner_(std::move(owner)),
    vis_key_(fieldname_),
    units_(),
    io_checkpoint_(true),
    io_vis_(true),
    initialized_(false)
{}


// pass-throughs for other functionality
void Record::WriteVis(const Visualization& vis) const {
  if (io_vis()) {
    data_.WriteVis(vis, vis_fieldname(), subfieldnames());
  }
}
void Record::WriteCheckpoint(const Checkpoint& chkp) const {
  if (io_checkpoint()) {
    data_.WriteCheckpoint(chkp, vis_fieldname());
  }
}
void Record::ReadCheckpoint(const Checkpoint& chkp) {
  if (io_checkpoint()) {
    data_.ReadCheckpoint(chkp, vis_fieldname());
  }
}

bool Record::Initialize(Teuchos::ParameterList& plist) {
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
    if (plist.isParameter("subfieldnames")) {
      set_subfieldnames(plist.get<Teuchos::Array<std::string> >("subfieldnames").toVector());
    }

    initialized_here = data_.Initialize(plist, vis_fieldname(), subfieldnames());
    if (initialized_here) set_initialized();
  }
  return initialized_here;
}


void Record::AssertOwnerOrDie(const Key& test_owner) const {
  if (test_owner != owner()) {
    Errors::Message message;
    message << "Record \"" << fieldname_ << "\" requested by \""
            << test_owner << "\" but is owned by \""
            << owner() << "\"";
    throw(message);
  }
};

} // namespace
