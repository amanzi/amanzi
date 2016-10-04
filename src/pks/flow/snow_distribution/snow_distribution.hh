/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#ifndef PK_FLOW_SNOW_DISTRIBUTION_HH_
#define PK_FLOW_SNOW_DISTRIBUTION_HH_

#include "upwinding.hh"

#include "Function.hh"

#include "Operator.hh"
#include "OperatorDiffusion.hh"
#include "OperatorAccumulation.hh"

#include "pk_factory.hh"
#include "pk_physical_bdf_base.hh"

namespace Amanzi {
namespace Flow {

class SnowDistribution : public PKPhysicalBDFBase {

public:
<<<<<<< HEAD
  SnowDistribution(Teuchos::Ptr<State> S, const Teuchos::RCP<Teuchos::ParameterList>& plist,
=======
  SnowDistribution(const Teuchos::RCP<Teuchos::ParameterList>& plist,
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
               Teuchos::ParameterList& FElist,
                   const Teuchos::RCP<TreeVector>& solution);
  
  // Virtual destructor
  virtual ~SnowDistribution() {}

  // main methods
  // -- Initialize owned (dependent) variables.
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {
    PKPhysicalBDFBase::commit_state(dt,S);
  }

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // error monitor
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du);

  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);
  
protected:
  // setup methods
  virtual void SetupSnowDistribution_(const Teuchos::Ptr<State>& S);
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

  // computational concerns in managing abs, rel perm
  // -- builds tensor K, along with faced-based Krel if needed by the rel-perm method
  virtual bool UpdatePermeabilityData_(const Teuchos::Ptr<State>& S);

  // boundary condition members
  virtual void UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S);
  
  // physical methods
  // -- diffusion term
  void ApplyDiffusion_(const Teuchos::Ptr<State>& S,const Teuchos::Ptr<CompositeVector>& g);
  // -- accumulation term
  void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g);

 protected:
  // control switches
  Operators::UpwindMethod upwind_method_;

  // coupling term
  bool full_jacobian_;

  // function for precip
  Teuchos::RCP<Function> precip_func_;
  
  // work data space
  Teuchos::RCP<Operators::Upwinding> upwinding_;

  // mathematical operators
  Teuchos::RCP<Operators::Operator> matrix_; // pc in PKPhysicalBDFBase
  Teuchos::RCP<Operators::OperatorDiffusion> matrix_diff_;
  Teuchos::RCP<Operators::OperatorDiffusion> preconditioner_diff_;
  Teuchos::RCP<Operators::OperatorAccumulation> preconditioner_acc_;

  // factory registration
  static RegisteredPKFactory<SnowDistribution> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
