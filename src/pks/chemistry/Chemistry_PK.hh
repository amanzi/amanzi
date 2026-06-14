/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*!

The chemistry header includes three parameters:

.. admonition:: chemistry_pk-spec

  * `"chemistry model`" ``[string]`` defines chemical model. The available options are `"Alquimia`"
    and `"Amanzi`" (default).

.. code-block:: xml

  <ParameterList name="_CHEMISTRY">
    <Parameter name="chemistry model" type="string" value="Amanzi"/>
  </ParameterList>


Geochemical engines
...................
Here we specify either the default or the third-party geochemical engine.

Common parameters
`````````````````
The following parameters are common for all supported engines.

.. admonition:: chemistry_params-spec

  * `"timestep control method`" ``[string]`` specifies timestep control method for chemistry subcycling.
    Choose either "fixed" (default) or "simple".  For option "fixed", timestep is fixed.
    For option "simple", the timestep is adjusted in response to stiffness of system of equations
    based on a simple scheme. This option require the following parameters: `"timestep cut threshold`",
    `"timestep cut factor`", `"timestep increase threshold`", and `"timestep increase factor`".

  * `"timestep cut threshold`" ``[int]`` is the number of Newton iterations that if exceeded
    will trigger a timestep cut. Default is 8.

  * `"max timestep (s)`" ``[double]`` is the maximum timestep that chemistry will allow the MPC to take.

  * `"initial timestep (s)`" ``[double]`` is the initial timestep that chemistry will ask the MPC to take.

  * `"timestep cut factor`" ``[double]`` is the factor by which the timestep is cut. Default is 2.0

  * `"timestep increase threshold`" ``[int]`` is the number of consecutive successful timesteps that
    will trigger a timestep increase. Default is 4.

  * `"timestep increase factor`" ``[double]`` is the factor by which the timestep is increased. Default is 1.2

  * `"free ion initial guess`" ``[double]`` provides an estimate of the free ion concentration for solutes.
    It used to help convergence of the initial solution of the chemistry. If this parameter is absent,
    a fraction (10%) of the total component concentration is used.

  * `"initial conditions time`" ``[double]`` specifies time for applying initial conditions. This parameter
    is useful for simulation restart. Default value is the state time when chemistry PK is instantiated.

*/

#ifndef AMANZI_CHEMISTRY_PK_HH_
#define AMANZI_CHEMISTRY_PK_HH_

#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "VerboseObject.hh"

// Amanzi
#ifdef ALQUIMIA_ENABLED
#include "ChemistryEngine.hh"
#endif
#include "Key.hh"
#include "Mesh.hh"
#include "PK_Physical_Default.hh"
#include "State.hh"

namespace Amanzi {

class TimestepController;

namespace AmanziChemistry {

class Chemistry_PK : public PK_Physical_Default {
 public:
  Chemistry_PK();
  Chemistry_PK(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& glist,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& soln);
  virtual ~Chemistry_PK() = default;

  // required members for PK interface
  virtual void parseParameterList() override;
  virtual void Initialize() override;

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

  virtual void set_dt(double dt) override {};
  virtual double get_dt() override;

  virtual void updateSubstate(const Tag& water_tag) = 0;

 protected:
  // error messages in Chemistry are rank-local.  This method does the
  // AllReduce to check errors across ranks.  Note that all three arguments are
  // accepted by reference, and updated with their equivalent GLOBAL values.
  void checkForError_(int& ierr, int& max_itrs, int& max_itrs_cell) const;

  // customization points for the beaker-level chemistry and substates
  virtual int advanceSingleCell_(int cell, double dt) = 0;
  virtual void copyFields_(const Tag& tag_dest, const Tag& tag_source) = 0;

 protected:
  int number_aqueous_components_;
  std::vector<std::string> aqueous_comp_names_;

  int number_gaseous_components_;
  std::vector<std::string> gaseous_comp_names_;

  int number_mineral_components_;
  std::vector<std::string> mineral_comp_names_;

  double saturation_tolerance_;

  // time controls
  double initial_conditions_time_;
  Teuchos::RCP<TimestepController> timestep_controller_;
  double dt_next_;
  int num_iterations_, num_successful_steps_;

  // discretization control
  Tag tcc_tag_current_;
  Tag tcc_tag_next_;
};

} // namespace AmanziChemistry
} // namespace Amanzi
#endif
