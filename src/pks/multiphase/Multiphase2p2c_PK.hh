/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Reimplementation of the original multiphase mode for testing.
  Solution vectors are pressure (l), saturation (l), and hydrogen 
  molar density (l).
*/

#ifndef AMANZI_MULTIPHASE_2P2C_PK_HH_
#define AMANZI_MULTIPHASE_2P2C_PK_HH_

// Amanzi
#include "Key.hh"
#include "PK_Factory.hh"

// Multiphase
#include "Multiphase_PK.hh"

namespace Amanzi {
namespace Multiphase {

class Multiphase2p2c_PK: public Multiphase_PK {
 public:
  Multiphase2p2c_PK(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& soln);

  ~Multiphase2p2c_PK() {};

  // modifying interface for PKs
  virtual void Setup(const Teuchos::Ptr<State>& S) override;
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) override;
  virtual std::string name() override { return "multiphase 2p2c"; }

  // interface multiphase models
  virtual void InitMPSolutionVector() override;
  virtual void InitMPPreconditioner() override;
  virtual void PopulateBCs(int icomp, bool flag) override;

  virtual std::pair<int, int> EquationToSolution(int neqn) override;
  virtual void ModifyEvaluators(int neqn) override;

 private:
  int missed_bc_faces_;

  Key advection_liquid_reduced_key_, molar_density_water_key_;
  Key ncp_f_key_, ncp_g_key_, ncp_fg_key_;

 private:
  // factory registration
  static RegisteredPKFactory<Multiphase2p2c_PK> reg_;
};

}  // namespace Multiphase
}  // namespace Amanzi
#endif
