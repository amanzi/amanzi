/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  MPC PK

  Weak coupling of mechanics and coupled flow PKs.
*/

#ifndef AMANZI_MECHANICS_COUPLED_FLOW_PK_HH_
#define AMANZI_MECHANICS_COUPLED_FLOW_PK_HH_

#include "Key.hh"

#include "PK_MPCWeak.hh"

namespace Amanzi {

class MechanicsFlowMatrixFracture_PK : public PK_MPCWeak {
 public:
  MechanicsFlowMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
                                 const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                 const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  virtual void Setup() override;

 private:
   Key hydrostatic_stress_key_;

 private:
  // factory registration
  static RegisteredPKFactory<MechanicsFlowMatrixFracture_PK> reg_;
};

} // namespace Amanzi

#endif
