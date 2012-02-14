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
#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "Field.hh"

namespace Amanzi {

class Field_CV : public Field {

public:
  // constructors
  Field_CV(std::string fieldname, std::string owner);
  Field_CV(std::string fieldname, std::string owner, Teuchos::RCP<CompositeVector>& data);

  // copy constructors
  explicit Field_CV(const Field_CV& other);
  Teuchos::RCP<Field> Clone() const;
  Teuchos::RCP<Field> Clone(std::string fieldname) const;
  Teuchos::RCP<Field> Clone(std::string fieldname, std::string owner) const;

  // destructor
  ~Field_CV() {}

  // subfield names
  std::vector< std::vector< std::string > > subfield_names() { return data_->subfield_names(); }
  void set_subfield_names(std::vector< std::vector< std::string > > subfield_names) {
    data_->set_subfield_names(subfield_names); }

  // data access and mutators
  Teuchos::RCP<const CompositeVector> GetFieldData() const { return data_; }
  Teuchos::RCP<CompositeVector> GetFieldData(std::string pk_name);

  void SetData(std::string pk_name, Teuchos::RCP<CompositeVector>& data);
  void SetData(std::string pk_name, const CompositeVector& data);

  // initialization
  virtual void Initialize(Teuchos::ParameterList& plist);

  // vis
  void WriteVis(Amanzi::Vis& vis);

protected:

  Teuchos::RCP<CompositeVector> data_;

private:
  // Assignment for the field disabled
  Field_CV& operator=(const Field_CV&);


}; // class Field

} // namespace Amanzi

#endif
