/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface of a Field containing a Scalar.

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#ifndef AMANZI_STATE_FIELD_SCALAR_HH_
#define AMANZI_STATE_FIELD_SCALAR_HH_

#include <string>
#include "Teuchos_RCP.hpp"

#include "Field.hh"

namespace Amanzi {

class Function;

class Field_Scalar : public Field {

public:
  // constructors
  Field_Scalar(std::string fieldname, std::string owner);
  Field_Scalar(std::string fieldname, std::string owner, Teuchos::RCP<double>& data);

  // copy constructor and assignment
  explicit Field_Scalar(const Field_Scalar& other);
  virtual Teuchos::RCP<Field> Clone() const;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname) const;
  virtual Teuchos::RCP<Field> Clone(std::string fieldname, std::string owner) const;

  // destructor
  ~Field_Scalar() {}

  // Creation of the actual data (data is created lazily, allowing empty fields).
  virtual void CreateData();

  // data access and mutators
  virtual Teuchos::RCP<const double> GetScalarData() const { return data_; }
  virtual Teuchos::RCP<double> GetScalarData();

  virtual void SetData(const Teuchos::RCP<double>& data);
  virtual void SetData(const double& data);

  // initialization
  virtual void Initialize(Teuchos::ParameterList& plist);

  // vis writing
  virtual void WriteVis(Visualization& vis);
  virtual void WriteCheckpoint(Checkpoint& ckp);

  // Compute from a function
  virtual void Compute(double time);

 protected:

  Teuchos::RCP<double> data_;
  Teuchos::RCP<Function> func_;

private:
  // operator= disabled
  Field_Scalar& operator=(const Field_Scalar&);

}; // class Field

} // namespace Amanzi

#endif
