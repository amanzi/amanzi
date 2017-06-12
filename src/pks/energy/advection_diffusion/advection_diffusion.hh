/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Start of a process kernel for the energy equation to be used in thermal
permafrost.  This starts with the simplification that T > T_freezing, limiting
us to the air-water system.
------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_DIFFUSION_HH_
#define PKS_ENERGY_DIFFUSION_HH_

#include "PK_Factory.hh"
#include "advection.hh"
#include "Operator.hh"
#include "OperatorDiffusion.hh"
#include "OperatorAdvection.hh"
#include "OperatorAccumulation.hh"
//#include "PK_PhysicalBDF_ATS.hh"
#include "pk_physical_bdf_default.hh"

namespace Amanzi {
namespace Energy {

class AdvectionDiffusion : public PK_PhysicalBDF_Default {

public:


  AdvectionDiffusion(Teuchos::ParameterList& FElist,
                     const Teuchos::RCP<Teuchos::ParameterList>& plist,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& solution) :
    PK(FElist, plist, S, solution),
    PK_PhysicalBDF_Default(FElist, plist, S, solution) {}


  // Virtual destructor
  virtual ~AdvectionDiffusion() {}

  // AdvectionDiffusion is a PK
  // -- Setup data.
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {
    PK_PhysicalBDF_Default::CommitStep(t_old, t_new, S);
  }

  // -- Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {}


  // AdvectionDiffusion is a BDFFnBase
  // computes the non-linear functional f = f(t,u,udot)
  virtual void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f);

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);


private:
  // helper methods for calling the above methods
  void AddAccumulation_(const Teuchos::RCP<CompositeVector> f);
  void AddAdvection_(const Teuchos::RCP<State> S,
                     const Teuchos::RCP<CompositeVector> f, bool negate);
  void ApplyDiffusion_(const Teuchos::RCP<State> S, const Teuchos::RCP<CompositeVector> f);

  // methods for applying/using the discretization/operators
  void UpdateBoundaryConditions_();
  void ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature);

  // misc setup information
  double dt_;
  bool implicit_advection_;

  // boundary conditions
  Teuchos::RCP<Functions::BoundaryFunction> bc_temperature_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;


  // mathematical operators
  Teuchos::RCP<Operators::OperatorDiffusion> matrix_diff_;
  Teuchos::RCP<Operators::OperatorAdvection> matrix_adv_;

  Teuchos::RCP<Operators::OperatorDiffusion> preconditioner_diff_;
  Teuchos::RCP<Operators::OperatorAccumulation> preconditioner_acc_;
  Teuchos::RCP<Operators::OperatorAdvection> preconditioner_adv_;

  // time integration
  double atol_;
  double rtol_;

  // factory registration
  static RegisteredPKFactory<AdvectionDiffusion> reg_;

};

} // namespace Energy
} // namespace Amanzi

#endif
