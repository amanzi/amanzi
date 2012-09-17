/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Start of a process kernel for the energy equation to be used in thermal
permafrost.  This starts with the simplification that T > T_freezing, limiting
us to the air-water system.
------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_TWOPHASE_HH_
#define PKS_ENERGY_TWOPHASE_HH_

#include "pk_factory.hh"
#include "bdf_time_integrator.hh"
#include "matrix_mfd.hh"

#include "pk_physical_bdf_base.hh"

namespace Amanzi {

// forward declarations
class MPCDiagonalFlowEnergy;
namespace Operators { class Advection; }
namespace Relations { class EOS; }
namespace Functions { class BoundaryFunction; }
namespace Energy { namespace EnergyRelations { class IEM; } }


namespace Energy {

class TwoPhase : public PKPhysicalBDFBase {

public:
  TwoPhase(Teuchos::ParameterList& plist, const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist,solution),
      PKPhysicalBDFBase(plist, solution) {
    plist_.set("solution key", "temperature");
  }

  // TwoPhase is a PK
  // -- Setup data
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);

  // -- Calculate any diagnostics prior to doing vis
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}


  // TwoPhase is a BDFFnBase
  // computes the non-linear functional f = f(t,u,udot)
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

protected:
  // for now, several points of entry into the science, as I'm not sure where
  // things will settle for a Phalanx-like system
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);
  virtual void SetupEnergy_(const Teuchos::Ptr<State>& S);

  virtual void UpdateBoundaryConditions_();
  virtual void ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature);
  virtual void DeriveFaceValuesFromCellValues_(const Teuchos::RCP<CompositeVector>& temp);

  virtual void AddAccumulation_(const Teuchos::RCP<CompositeVector> f);
  virtual void AddAdvection_(const Teuchos::RCP<State> S,
                     const Teuchos::RCP<CompositeVector> f, bool negate);
  virtual void ApplyDiffusion_(const Teuchos::RCP<State> S,
          const Teuchos::RCP<CompositeVector> f);

protected:
  // boundary conditions
  Teuchos::RCP<Functions::BoundaryFunction> bc_temperature_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;
  std::vector<Operators::Matrix_bc> bc_markers_;
  std::vector<double> bc_values_;

  // operators
  Teuchos::RCP<Operators::Advection> advection_;
  Teuchos::RCP<Operators::MatrixMFD> matrix_;
  Teuchos::RCP<Operators::MatrixMFD> preconditioner_;

  // models for evaluating BCs
  Teuchos::RCP<Relations::EOS> eos_liquid_;
  Teuchos::RCP<EnergyRelations::IEM> iem_liquid_;

private:
  // factory registration
  static RegisteredPKFactory<TwoPhase> reg_;

  // Energy has a friend in couplers...
  friend class Amanzi::MPCDiagonalFlowEnergy;

};

} // namespace Energy
} // namespace Amanzi

#endif
