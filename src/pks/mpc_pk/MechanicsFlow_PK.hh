/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  MPC PK

  Iterative coupling of mechanics and flow PKs via undrained split.
*/

#ifndef AMANZI_MECHANICS_FLOW_PK_HH_
#define AMANZI_MECHANICS_FLOW_PK_HH_

#include "Key.hh"
#include "PorosityEvaluator.hh"

#include "PK_MPCSequential.hh"

namespace Amanzi {

// class MechanicsFlow_PK : public PK_MPCWeak {
class MechanicsFlow_PK : public PK_MPCSequential {
 public:
  MechanicsFlow_PK(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  virtual void Setup() override;
  virtual void Initialize() override;
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;

 private:
  void EvaluateForDarcy_();
  void EvaluateForRichards_();

 private:
  Key domain_;
  const Teuchos::RCP<Teuchos::ParameterList> glist_;

  Key displacement_key_, hydrostatic_stress_key_, vol_strain_key_;
  Key pressure_key_, saturation_liquid_key_, water_storage_key_, porosity_key_;
  Key undrained_split_coef_key_;

  bool thermal_flow_;

  static RegisteredPKFactory<MechanicsFlow_PK> reg_;
};

} // namespace Amanzi

#endif
