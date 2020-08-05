/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface of a Field containing a vector that is constant over the mesh.

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
  Field_ConstantVector(std::string fieldname, std::string owner, int dimension);
  Field_ConstantVector(std::string fieldname, std::string owner,
                       Teuchos::RCP<Epetra_Vector>& data);

  // copy constructor and assignment
  explicit Field_ConstantVector(const Field_ConstantVector& other);
  virtual Teuchos::RCP<Field> Clone() const;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname) const;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname, std::string owner) const;

  // destructor
  ~Field_ConstantVector() {}

  // Creation of the actual data (data is created lazily, allowing empty fields).
  virtual void CreateData();

  // dimension
  int dimension() { return dimension_; }
  void set_dimension(int dimension);

  // names for the vector components
  std::vector<std::string> subfield_names() { return subfield_names_; }
  void set_subfield_names(std::vector<std::string> subfield_names) {
    subfield_names_ = subfield_names; }

  // data access and mutators
  Teuchos::RCP<const Epetra_Vector> GetConstantVectorData() const { return data_; }
  Teuchos::RCP<Epetra_Vector> GetConstantVectorData();

  void SetData(const Teuchos::RCP<Epetra_Vector>& data);
  void SetData(const Epetra_Vector& data);

  // initialization
  virtual void Initialize(Teuchos::ParameterList& plist);

  // visualization
  void WriteVis(Visualization& vis);

  // checkpoint
  void WriteCheckpoint(Checkpoint& ckp);

protected:

  Teuchos::RCP<Epetra_Vector> data_;
  std::vector<std::string> subfield_names_;
  int dimension_;

private:
  // operator= disabled
  Field_ConstantVector& operator=(const Field_ConstantVector&);

}; // class Field

} // namespace Amanzi

#endif
