/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for a Field containing a CompositeVector.  Field is not intended so
much to hide implementation of data as to restrict write access to data.  It
freely passes out pointers to its private data, but only passes out read-only
const pointers unless you have the secret password (AKA the name of the
process kernel that owns the data).

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
#include "Field.hh"

namespace Amanzi {

class Field_CV : public Field {

public:
  // constructors
  Field_CV(std::string fieldname, std::string owner="state");

  // copy constructor and assignment
  Teuchos::RCP<Field_CV> Clone() const;
  operator=(const Field&);

  // data access and mutators
  virtual Teuchos::RCP<const CompositeVector> data() const { return data_; }
  virtual Teuchos::RCP<CompositeVector> data(std::string pk_name) {
    assert_owner_or_die_(pk_name);
    return data_;
  }

  virtual void set_data(std::string pk_name, Teuchos::RCP<CompositeVector>& data) {
    assert_owner_or_die_(pk_name);
    data_ = data;
  }

private:

  Teuchos::RCP<CompositeVector> data_;

}; // class Field

} // namespace Amanzi

#endif
