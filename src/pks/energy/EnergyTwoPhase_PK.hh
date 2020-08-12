/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Process kernel for the thermal Richards flow.
*/

#ifndef AMANZI_ENERGY_TWOPHASE_PK_HH_
#define AMANZI_ENERGY_TWOPHASE_PK_HH_

// Amanzi
#include "BDF1_TI.hh"
#include "EOS.hh"
#include "IEM.hh"
#include "PK_Factory.hh"

// Energy
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

  // methods required for PK intrefcae
  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false);
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {};

  double get_dt() { return dt_; }
  void set_dt(double dt) { dt_ = dt; }

  virtual std::string name() { return "two-phase energy"; }

  // methods required for time integration interface
  // -- computes the non-linear functional f = f(t,u,udot) and related norm.
  virtual void FunctionalResidual(const double t_old, double t_new,
                          Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new,
                          Teuchos::RCP<TreeVector> g);
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);

  // -- management of the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double dt);

  // access method for unit tests
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > bdf1_dae() { return bdf1_dae_; }

 private:
  void InitializeFields_();

 protected:
  // models for evaluating total energy
  Teuchos::RCP<AmanziEOS::EOS> eos_liquid_;
  Teuchos::RCP<IEM> iem_liquid_;

 private:
  // primary field
  const Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<CompositeVector> solution;

  // time stepping
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > bdf1_dae_;
  int num_itrs_;
  double dt_, dt_next_;

  // factory registration
  static RegisteredPKFactory<EnergyTwoPhase_PK> reg_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
