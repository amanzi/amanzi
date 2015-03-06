/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Process kernel for thermal Richards' flow.
*/

#ifndef AMANZI_ENERGY_TWOPHASE_PK_HH_
#define AMANZI_ENERGY_TWOPHASE_PK_HH_

#include "eos.hh"
#include "iem.hh"
#include "PK_Factory.hh"
#include "Energy_PK.hh"

namespace Amanzi {
namespace Energy {

class EnergyTwoPhase_PK : public Energy_PK {

public:
  EnergyTwoPhase_PK(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& glist,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& soln);
  virtual ~EnergyTwoPhase_PK() {};

  // Required PK members.
  virtual void Setup();
  virtual void Initialize();
  virtual std::string name() { return "two-phase energy"; }
  virtual void CommitStep(double t_old, double t_new);

  virtual void Functional(const double t_old, double t_new,
                          Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new,
                          Teuchos::RCP<TreeVector> g);
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double dt);

  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);

 protected:
  // models for evaluating enthalpy
  Teuchos::RCP<Relations::EOS> eos_liquid_;
  Teuchos::RCP<IEM> iem_liquid_;

 private:
  Teuchos::RCP<Teuchos::ParameterList> ep_list_;

  // primary field
  const Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<CompositeVector> solution;

private:
  // factory registration
  static RegisteredPKFactory<EnergyTwoPhase_PK> reg_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
