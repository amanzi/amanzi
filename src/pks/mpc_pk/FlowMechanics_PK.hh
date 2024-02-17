/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  MPC PK

  Sequential coupling of flow and mechanics PKs via the fixed stress 
  split algorithm.
*/

#ifndef AMANZI_FLOW_MECHANICS_PK_HH_
#define AMANZI_FLOW_MECHANICS_PK_HH_

#include "Key.hh"

#include "PK_MPCSequential.hh"

namespace Amanzi {

class FlowMechanics_PK : public PK_MPCSequential {
 public:
  FlowMechanics_PK(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  virtual void Setup() override;
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;
  virtual double
  ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) override;
  virtual void CommitSequentialStep(Teuchos::RCP<const TreeVector> u_old,
                                    Teuchos::RCP<const TreeVector> u_new) override;

 private:
  Key domain_;
  const Teuchos::RCP<Teuchos::ParameterList> glist_;

  Key displacement_key_, hydrostatic_stress_key_, vol_strain_key_;
  Key pressure_key_, porosity_key_, saturation_liquid_key_, water_storage_key_;

  bool thermal_flow_;

 private:
  static RegisteredPKFactory<FlowMechanics_PK> reg_;
};

} // namespace Amanzi

#endif
