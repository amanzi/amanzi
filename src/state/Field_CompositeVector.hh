/*
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
// Interface for a Field containing a CompositeVector.

/*

  Field_CompositeVector is metadata for a CompositeVector, but also includes
functionality to initialize, restart, visualize, etc the data.  Parameters
available here are provided for use in InitialConditions_ specs.

*/

/*!

The vast majority of Fields are vectors of data on a mesh entity.  This is a
CompositeVector in Amanzi-ATS (as it may consist of multiple components, each
on its own mesh entity).  Initializing these components may be done in avariety
of ways:


``[initial-conditions-compositevector-spec]``

* `"restart file`" ``[string]`` **optional** If provided, read IC value from a
  checkpoint file of this name.

* `"cells from file`" ``[string]`` **optional** If provided, read IC value for
  cells from a checkpoint file of this name.

* `"exodus file initialization`" ``[exodus-file-reader-spec]`` **optional** If provided, initialize data from an Exodus file.


* `"constant`" ``[double]`` **optional** If provided, set the IC to this value everywhere.

* `"initialize from 1D column`" ``[initialize-1d-column-spec]`` **optional**
   If provided, initialize data by interpolating from a 1D solution in depth
   relative to the surface.

* `"function`" ``[composite-vector-function-spec-list]`` **optional** If
   provided, a region-based list of functions to evaluate piecewise over the
   domain.

* `"initialize faces from cell`" ``[bool]`` **false** In the case of mimetic
   or other face-inclusive discretizations, calculate initial face data by
   averages of the neighboring cells.

* `"write checkpoint`" ``[bool]`` **true** Write this data when checkpointing.

* `"write vis`" ``[bool]`` **true** Write this data when visualizing.


Capability for 1D column initialization:

``[initialize-1d-column-spec]``
* `"file`" ``[string]`` Filename including data.
* `"z header`" ``[string]`` **"/z"**  Name of the data in the above file for z-coordinate of this field.
* `"f header`" ``[string]`` **"/field-name"** Name of the data in the above file for this field.
* `"coordinate orientation`" ``[string]`` **"standard"** Either `"standard`" or `"depth`" -- is positive z up (standard) or down (depth)?

ONE OF:
* `"surface sideset`" ``[string]`` Name of the face-set in the mesh which determines the zero surface for interpolation.  This must cover the entire top surface boundary.
OR:
* `"surface sidesets`" ``[Array(string)]`` List of the names for the case of multiple sidesets needed.
END

ExodusII file reader control:

``[exodus-file-reader-spec]``
* `"file`" ``[string]`` Filename
* `"attributes`" ``[Array(string)]`` List of Exodus attributes to read, unclear when/if/how this works?

*/

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
  void WriteVis(Visualization& vis);
  void WriteCheckpoint(Checkpoint& ckp);
  bool ReadCheckpoint(HDF5_MPI& file);

  long int GetLocalElementCount();

  // initialization
  // -- this function assumes that all components of the data vector
  //    were initialized. If component is missing, it returns false
  virtual bool initialized() const;
  // -- sets all components of the data vector to the given flag
  virtual void set_initialized(bool initialized=true);

  // component-wise initialization
  // -- returns false if component is either missing or not initialized
  bool initialized(const std::string& comp) const;
  // -- does nothing if component does not exist in the data vector
  void set_initialized(const std::string& comp, bool initialized=true);

protected:
  bool ReadCheckpoint_(std::string filename);
  void ReadCellsFromCheckpoint_(const std::string& filename); // for ICs
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

  // component information
  std::map<std::string, bool> initialized_comp_;
};

} // namespace Amanzi

#endif
