/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for a Field.  Field is not intended so much to hide implementation
of data as to restrict write access to data.  It freely passes out pointers to
its private data, but only passes out read-only const pointers unless you have
the secret password (AKA the name of the process kernel that owns the data).

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#ifndef STATE_FIELD_CV_HH_
#define STATE_FIELD_CV_HH_

#include <string>
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

#include "Mesh.hh"
#include "CompositeVector.hh"

namespace Amanzi {

class Field {

public:
  // default constructor, does not generate data
  Field(std::string fieldname, std::string owner="state");
  Field(std::string fieldname, std::string owner="state",
        Teuchos::RCP<CompositeVector>& data);

  // copy constructor and assignment
  Field(const Field&);
  Field& operator=(const Field&);

  // access
  std::string fieldname() const { return fieldname_; }
  std::string owner() const { return owner_; }
  std::vector<std::vector< std::string > > subfield_names() const { return subfield_names_};
  bool io_restart() const { return io_restart_; }
  bool io_vis() const { return io_vis_; }
  bool initialized() const { return initialized_; }

  Teuchos::RCP<const CompositeVector> data() const { return data_; }
  Teuchos::RCP<CompositeVector> data(std::string pk_name);

  // mutators
  void set_owner(std::string owner) { owner_ = owner; }
  void set_io_restart(bool io_restart=true) { io_restart_ = io_restart; }
  void set_io_vis(bool io_vis=true) { io_vis_ = io_vis; }
  void set_initialized(bool initialized=true) { initialized_ = initialized; }

  void set_subfield_names(std::vector< std::vector<std::string> > subfieldnames) {
    subfieldnames_ = subfieldnames;
  }
  void set_subfield_names(std::vector<std::string> subfieldnames) {
    subfieldnames_.PushBack(subfieldnames);
  }

  void set_data(std::string pk_name, Teuchos::RCP<CompositeVector>& data);

private:

  // consistency checking
  void assert_owner_or_die_(std::string pk_name) const;

  Teuchos::RCP<CompositeVector> data_;
  std::string fieldname_;
  std::string owner_;
  std::vector< std::vector<std::string> > subfieldnames_;

  bool io_restart_;
  bool io_vis_;
  bool initialized_;

}; // class Field

} // namespace Amanzi

#endif
