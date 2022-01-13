/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
           Konstantin Lipnikov

  Process kernel that strongly couples Flow PK with Energy PK.
*/

#ifndef AMANZI_FLOW_ENERGY_PK_HH_
#define AMANZI_FLOW_ENERGY_PK_HH_

#include "Teuchos_RCP.hpp"

#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "PK_BDF.hh"
#include "PK_MPCStrong.hh"
#include "PK_Factory.hh"

namespace Amanzi {

class FlowEnergy_PK : public PK_MPCStrong<PK_BDF> {
 public:
  FlowEnergy_PK(Teuchos::ParameterList& pk_tree,
                const Teuchos::RCP<Teuchos::ParameterList>& glist,
                const Teuchos::RCP<State>& S,
                const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  virtual void Setup() override;
  virtual void Initialize() override;  

  // -- dt is the minimum of the sub pks
  // virtual double get_dt();
  // virtual void set_dt(double dt);

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;

  std::string name() override { return "thermal flow"; }

 private:
  const Teuchos::RCP<Teuchos::ParameterList>& glist_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Key domain_;  // computational domain

  Teuchos::RCP<EvaluatorIndependentFunction> particle_density_eval;
  Teuchos::RCP<EvaluatorIndependentFunction> porosity_eval;
  Teuchos::RCP<EvaluatorIndependentFunction> saturation_liquid_eval;

  // keys
  Key ie_rock_key_, ie_gas_key_, ie_liquid_key_, energy_key_, prev_energy_key_;
  Key particle_density_key_;
  Key mol_density_liquid_key_, mol_density_gas_key_, mass_density_liquid_key_;
  Key pressure_key_, sat_liquid_key_, prev_sat_liquid_key_;
  Key wc_key_, prev_wc_key_;
  Key viscosity_liquid_key_;

  // eos
  std::string eos_table_;

  // factory registration
  static RegisteredPKFactory<FlowEnergy_PK> reg_;
};

}  // namespace Amanzi
#endif
