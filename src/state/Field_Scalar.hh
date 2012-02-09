/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface of a Field containing a Scalar.  Field is not intended so
much to hide implementation of data as to restrict write access to data.  It
freely passes out pointers to its private data, but only passes out read-only
const pointers unless you have the secret password (AKA the name of the
process kernel that owns the data).

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#ifndef STATE_FIELD_SCALAR_HH_
#define STATE_FIELD_SCALAR_HH_

#include <string>
#include "Teuchos_RCP.hpp"

#include "Field.hh"

namespace Amanzi {

class Field_Scalar : public Field {

public:
  // constructors
  Field_Scalar(std::string fieldname);
  Field_Scalar(std::string fieldname, std::string owner);
  Field_Scalar(std::string fieldname, Teuchos::RCP<double>& data);
  Field_Scalar(std::string fieldname, std::string owner, Teuchos::RCP<double>& data);

  // copy constructor and assignment
  explicit Field_Scalar(const Field_Scalar& other);
  virtual Teuchos::RCP<Field> Clone() const;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname) const;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname, std::string owner) const;

  // data access and mutators
  Teuchos::RCP<const double> scalar_data() const { return data_; }
  Teuchos::RCP<double> scalar_data(std::string pk_name);

  void set_data(std::string pk_name, Teuchos::RCP<double>& data);

protected:

  Teuchos::RCP<double> data_;

private:
  // operator= disabled
  Field_Scalar& operator=(const Field_Scalar&);

}; // class Field

} // namespace Amanzi

#endif
