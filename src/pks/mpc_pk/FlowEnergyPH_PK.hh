/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
           Konstantin Lipnikov
*/

/*!

Process kernel that strongly couples Flow PK with Energy PK.

*/

#ifndef AMANZI_FLOW_ENERGY_PRESSURE_ENTHALPY_PK_HH_
#define AMANZI_FLOW_ENERGY_PRESSURE_ENTHALPY_PK_HH_

#include "Teuchos_RCP.hpp"

#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "PDE_Accumulation.hh"
#include "PDE_Advection.hh"
#include "PDE_Diffusion.hh"
#include "PK_BDF.hh"
#include "PK_MPCStrong.hh"
#include "PK_Factory.hh"
#include "ThermodynamicStateEvaluators.hh"

namespace Amanzi {

class FlowEnergyPH_PK : public PK_MPCStrong<PK_BDF> {
 public:
  FlowEnergyPH_PK(Teuchos::ParameterList& pk_tree,
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

  // -- computes non-linear functional f = f(t,u)
  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<const TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> f) override;

  // -- preconditioner
  virtual void UpdatePreconditioner(double t,
                                    Teuchos::RCP<const TreeVector> up,
                                    double dt) override;

  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> Pu) override;

  // -- error norm for coupled system
  std::string name() override { return "flow and energy ph"; }

  // -- L-scheme for flow equaiton
  virtual std::vector<Key> SetupLSchemeKey(Teuchos::ParameterList& plist) override;

 private:
  const Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Key domain_; // computational domain

  Teuchos::RCP<Operators::Operator> op10_, op01_;
  Teuchos::RCP<Operators::PDE_Diffusion> pde10_diff_;
  Teuchos::RCP<Operators::PDE_Advection> pde10_adv_;
  Teuchos::RCP<Operators::PDE_Accumulation> pde10_acc_, pde01_acc_;
  bool symbolic_assembly_complete_ = false;

  // keys
  Key pressure_key_, enthalpy_key_, temperature_key_;
  Key energy_key_, particle_density_key_, conductivity_key_;
  Key state_key_, viscosity_liquid_key_, mol_density_liquid_key_, iso_compressibility_key_;
  Key mol_flowrate_key_, water_storage_key_;
  Key  bcs_flow_key_, bcs_enthalpy_key_;

  // factory registration
  static RegisteredPKFactory<FlowEnergyPH_PK> reg_;
};

} // namespace Amanzi
#endif
