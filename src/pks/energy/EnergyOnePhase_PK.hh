/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Process kernel for single phase energy equation.
*/

#ifndef AMANZI_ENERGY_ONEPHASE_PK_HH_
#define AMANZI_ENERGY_ONEPHASE_PK_HH_

// Amanzi
#include "BDF1_TI.hh"
#include "IEM.hh"
#include "PK_Factory.hh"

// Energy
#include "Energy_PK.hh"

namespace Amanzi {
namespace Energy {

class EnergyOnePhase_PK : public Energy_PK {
 public:
  EnergyOnePhase_PK(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& glist,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& soln);
  virtual ~EnergyOnePhase_PK(){};

  // methods required for PK intrefcae
  virtual void Setup() final;
  virtual void Initialize() final;

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) final;
  virtual void CalculateDiagnostics(const Tag& tag) final{};

  double get_dt() final { return dt_; }
  void set_dt(double dt) final { dt_ = dt; }

  virtual std::string name() { return "one-phase energy"; }

  // methods required for time integration interface
  // -- computes the non-linear functional f = f(t,u,udot) and related norm.
  virtual void FunctionalResidual(const double t_old,
                                  double t_new,
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> g) final;
  virtual double
  ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) final;

  // -- management of the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double dt);

  // access method for unit tests
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace>> bdf1_dae() { return bdf1_dae_; }

 private:
  void InitializeFields_();

 private:
  // primary field
  const Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<CompositeVector> solution;

  // time stepping
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace>> bdf1_dae_;
  int num_itrs_;
  double dt_, dt_next_;

  // factory registration
  static RegisteredPKFactory<EnergyOnePhase_PK> reg_;
};

} // namespace Energy
} // namespace Amanzi

#endif
