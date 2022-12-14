/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*
  PK for coupling of Flow PK with Transport_PK and Chemestry_PK

*/


#ifndef AMANZI_PRESSURESATURATION_PK_HH_
#define AMANZI_PRESSURESATURATION_PK_HH_

#include "Teuchos_RCP.hpp"

#include "PK.hh"
#include "PK_Factory.hh"
#include "MPCSubcycled.hh"

namespace Amanzi {

class PressureSaturation_PK : public MPCSubcycled {
 public:
  PressureSaturation_PK(Teuchos::ParameterList& pk_tree,
                        const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                        const Teuchos::RCP<State>& S,
                        const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt);

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

  // virtual void Initialize();

  std::string name() { return "pressure saturation"; }

 private:
  // factory registration
  static RegisteredPKFactory<PressureSaturation_PK> reg_;
};

} // namespace Amanzi
#endif
