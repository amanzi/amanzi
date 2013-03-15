/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for a Field containing a CompositeVector.

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#ifndef STATE_FIELD_CV_HH_
#define STATE_FIELD_CV_HH_

#include <string>
#include "Teuchos_RCP.hpp"

#include "composite_vector.hh"
#include "field.hh"

namespace Amanzi {

class FieldCompositeVector : public Field {

public:
  // constructors
  FieldCompositeVector(std::string fieldname, std::string owner);
  FieldCompositeVector(std::string fieldname, std::string owner, Teuchos::RCP<CompositeVector>& data);

  // copy constructors
  explicit FieldCompositeVector(const FieldCompositeVector& other);
  Teuchos::RCP<Field> Clone() const;
  Teuchos::RCP<Field> Clone(std::string fieldname) const;
  Teuchos::RCP<Field> Clone(std::string fieldname, std::string owner) const;

  // destructor
  ~FieldCompositeVector() {}

  // data creation
  void CreateData();

  // subfield names
  std::vector<std::vector<std::string > > subfield_names() { return subfield_names_; }
  void set_subfield_names(std::vector<std::vector<std::string> >& subfield_names) {
    subfield_names_ = subfield_names; }

  // data access and mutators
  Teuchos::RCP<const CompositeVector> GetFieldData() const { return data_; }
  Teuchos::RCP<CompositeVector> GetFieldData();

  void SetData(const Teuchos::RCP<CompositeVector>& data);
  void SetData(const CompositeVector& data);

  // initialization
  virtual void Initialize(Teuchos::ParameterList& plist);

  // visualization and checkpoint i/o
  void WriteVis(const Teuchos::Ptr<Visualization>& vis);
  void WriteCheckpoint(const Teuchos::Ptr<Checkpoint>& ckp);
  void ReadCheckpoint(const Teuchos::Ptr<HDF5_MPI>& file);

protected:
  void ReadCheckpoint_(std::string filename);
  void ReadCellsFromCheckpoint_(std::string filename); // for ICs

  Teuchos::RCP<CompositeVector> data_;
  std::vector<std::vector<std::string> > subfield_names_;

private:
  // Assignment for the field disabled
  FieldCompositeVector& operator=(const FieldCompositeVector&);

  // check to ensure subfield names are set
  void EnsureSubfieldNames_();

}; // class Field

} // namespace Amanzi

#endif
