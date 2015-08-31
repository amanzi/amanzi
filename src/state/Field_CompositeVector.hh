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

#include "CompositeVector.hh"
#include "Field.hh"

namespace Amanzi {

void DeriveFaceValuesFromCellValues(CompositeVector& cv);


class Field_CompositeVector : public Field {

public:
  // constructors
  Field_CompositeVector(std::string fieldname, std::string owner);
  Field_CompositeVector(std::string fieldname, std::string owner,
                        const std::vector<std::vector<std::string> >& subfield_names);
  Field_CompositeVector(std::string fieldname, std::string owner, Teuchos::RCP<CompositeVector>& data);

  // copy constructors
  explicit Field_CompositeVector(const Field_CompositeVector& other);
  Teuchos::RCP<Field> Clone() const;
  Teuchos::RCP<Field> Clone(std::string fieldname) const;
  Teuchos::RCP<Field> Clone(std::string fieldname, std::string owner) const;

  // destructor
  ~Field_CompositeVector() {}

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
  bool ReadCheckpoint(const Teuchos::Ptr<HDF5_MPI>& file);

  long int GetLocalElementCount();

protected:
  bool ReadCheckpoint_(std::string filename);
  void ReadCellsFromCheckpoint_(std::string filename); // for ICs
  void ReadAttributeFromExodusII_(Teuchos::ParameterList& plist); 
  void ReadVariableFromExodusII_(Teuchos::ParameterList& plist); 
  void InitializeFromColumn_(Teuchos::ParameterList& plist);

  Teuchos::RCP<CompositeVector> data_;
  std::vector<std::vector<std::string> > subfield_names_;

private:
  // Assignment for the field disabled
  Field_CompositeVector& operator=(const Field_CompositeVector&);

  // check to ensure subfield names are set
  void EnsureSubfieldNames_();

}; // class Field

} // namespace Amanzi

#endif
