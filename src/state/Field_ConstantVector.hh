/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface of a Field containing a vector that is constant over the mesh.
Field is not intended so much to hide implementation of data as to restrict
write access to data.  It freely passes out pointers to its private data, but
only passes out read-only const pointers unless you have the secret password
(AKA the name of the process kernel that owns the data).

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#ifndef STATE_FIELD_CONSTANTVECTOR_HH_
#define STATE_FIELD_CONSTANTVECTOR_HH_

#include <string>
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

#include "Field.hh"

namespace Amanzi {

class Field_ConstantVector : public Field {

public:
  // constructors
  Field_ConstantVector(std::string fieldname, std::string owner);
  Field_ConstantVector(std::string fieldname, std::string owner,
                       Teuchos::RCP<Epetra_Vector>& data);

  // copy constructor and assignment
  explicit Field_ConstantVector(const Field_ConstantVector& other);
  virtual Teuchos::RCP<Field> Clone() const;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname) const;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname, std::string owner) const;

  // data access and mutators
  Teuchos::RCP<const Epetra_Vector> constant_vector_data() const { return data_; }
  Teuchos::RCP<Epetra_Vector> constant_vector_data(std::string pk_name);

  void set_data(std::string pk_name, Teuchos::RCP<Epetra_Vector>& data);
  void set_data(std::string pk_name, const Epetra_Vector& data);

  // initialization
  virtual void Initialize(Teuchos::ParameterList& plist);

protected:

  Teuchos::RCP<Epetra_Vector> data_;

private:
  // operator= disabled
  Field_ConstantVector& operator=(const Field_ConstantVector&);

}; // class Field

} // namespace Amanzi

#endif
