/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

The Amanzi chemistry process kernel uses the following parameters.

* `"thermodynamic database`" [list]

  * `"file`" [string] is the name of the chemistry database file, relative to the execution directory.

  * `"format`" [string] is the format of the database file. Actual database format is not XML and
    is the same as described for the 2010 demo with additions for the new chemical processes.
    Valid values: "simple".

* `"minerals`" [Array(string)] is the list of mineral names.

* `"sorption sites`" [Array(string)] is the list of sorption sites.

* `"activity model`" [string] is the type of model used for activity corrections.
  Valid options are `"unit`", `"debye-huckel`", and `"pitzer-hwm`",

* `"tolerance`" [double] defines tolerance in Newton solves inside the chemistry library.

* `"maximum Newton iterations`" [int] is the maximum number of iteration the chemistry
  library can take.

* `"auxiliary data`" [Array(string)] defines additional chemistry related data that the user
  can request be saved to vis files. Currently `"pH`" is the only variable supported.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="_CHEMISTRY">
    <ParameterList name="thermodynamic database">
      <Parameter name="file" type="string" value="_TRITIUM.bgd"/>
      <Parameter name="format" type="string" value="simple"/>
    </ParameterList>
    <Parameter name="activity model" type="string" value="unit"/>
    <Parameter name="tolerance" type="double" value="1.5e-12"/>
    <Parameter name="maximum Newton iterations" type="int" value="25"/>
    <Parameter name="max time step (s)" type="double" value="1.5e+07"/>
    <Parameter name="auxiliary data" type="Array(string)" value="{pH}"/>
    <Parameter name="number of component concentrations" type="int" value="1"/>
    <Parameter name="time step control method" type="string" value="simple"/>
  </ParameterList>
  </ParameterList>

*/

#ifndef CHEMISTRY_AMANZI_PK_HH_
#define CHEMISTRY_AMANZI_PK_HH_

#include <string>
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "Beaker.hh"
#include "BeakerFields.hh"
#include "BeakerState.hh"
#include "Chemistry_PK.hh"
#include "Key.hh"
#include "PK_Factory.hh"
#include "Mesh.hh"
#include "TreeVector.hh"

namespace Amanzi {
namespace AmanziChemistry {

// Trilinos based chemistry process kernel for the unstructured mesh
class Amanzi_PK : public Chemistry_PK {
 public:
  Amanzi_PK(Teuchos::ParameterList& pk_tree,
            const Teuchos::RCP<Teuchos::ParameterList>& glist,
            const Teuchos::RCP<State>& S,
            const Teuchos::RCP<TreeVector>& soln);

  // members required by PK interface
  virtual void Setup() final;
  virtual void Initialize() final;

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) final;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) final;
  virtual void CalculateDiagnostics(const Tag& tag) final { extra_chemistry_output_data(); }

  virtual std::string name() { return "chemistry amanzi"; }

  // The following two routines provide the interface for
  // output of auxillary cellwise data from chemistry
  Teuchos::RCP<Epetra_MultiVector> extra_chemistry_output_data();
  void set_chemistry_output_names(std::vector<std::string>* names);

  // functions used in Rransport PK
  void CopyCellStateToBeakerState(int c, Teuchos::RCP<Epetra_MultiVector> aqueous_components);

  // access
  std::shared_ptr<Beaker> get_engine() { return chem_; }
  const BeakerParameters& beaker_parameters() const { return beaker_parameters_; }
  BeakerState beaker_state() { return beaker_state_; }

 private:
  void AllocateAdditionalChemistryStorage_();

  void XMLParameters();
  void SetupAuxiliaryOutput();

  void InitializeBeakerFields_();

  void CopyBeakerStructuresToCellState(int c, Teuchos::RCP<Epetra_MultiVector> aqueous_components);

  void EstimateNextTimeStep_(double t_old, double t_new);

 protected:
  Teuchos::RCP<TreeVector> soln_;

 private:
  std::shared_ptr<Beaker> chem_;
  BeakerParameters beaker_parameters_;
  BeakerState beaker_state_, beaker_state_copy_;
  BeakerFields bf_;

  Key prev_saturation_key_; // move to base class ???

  std::string dt_control_method_;
  double dt_int_, dt_global_; // interpolation and global times

  std::vector<std::string> aux_names_;
  std::vector<int> aux_index_;
  Teuchos::RCP<Epetra_MultiVector> aux_data_;

  int ncells_owned_;

 private:
  // factory registration
  static RegisteredPKFactory<Amanzi_PK> reg_;
};

} // namespace AmanziChemistry
} // namespace Amanzi
#endif
