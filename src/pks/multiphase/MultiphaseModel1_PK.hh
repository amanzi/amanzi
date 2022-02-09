/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Solution vector: pressure (l), saturation (l), mole fraction (g).
*/

#ifndef AMANZI_MULTIPHASE_MODEL1_PK_HH_
#define AMANZI_MULTIPHASE_MODEL1_PK_HH_

// Amanzi
#include "Key.hh"
#include "PK_Factory.hh"

// Multiphase
#include "Multiphase_PK.hh"

namespace Amanzi {
namespace Multiphase {

class MultiphaseModel1_PK: public Multiphase_PK {
 public:
  MultiphaseModel1_PK(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln);

  ~MultiphaseModel1_PK() {};

  // modifying interface for PKs
  virtual void Setup() override;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;
  virtual std::string name() override { return "multiphase pl sl xg"; }

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du) override; 

  // interface multiphase models
  virtual void PopulateBCs(int icomp, bool flag) override;

  virtual SolutionStructure EquationToSolution(int neqn) override;
  virtual void ModifyEvaluators(int neqn) override;

 private:
  int missed_bc_faces_;

  Key advection_water_key_, pressure_vapor_key_, x_vapor_key_;
  Key diffusion_liquid_key_, diffusion_gas_key_, diffusion_vapor_key_; 
  Key molecular_diff_liquid_key_, molecular_diff_gas_key_; 

 private:
  // factory registration
  static RegisteredPKFactory<MultiphaseModel1_PK> reg_;
};

}  // namespace Multiphase
}  // namespace Amanzi
#endif
