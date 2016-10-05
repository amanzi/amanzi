/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Base energy PK.

This provides the base of an advection-diffusion equation for energy.

------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_BASE_HH_
#define PKS_ENERGY_BASE_HH_

#include "pk_factory.hh"

#include "OperatorDiffusion.hh"
#include "OperatorDiffusionMFD.hh"
#include "OperatorAccumulation.hh"
#include "OperatorAdvection.hh"

#include "pk_physical_bdf_base.hh"
#include "upwinding.hh"

namespace Amanzi {

// forward declarations
namespace Operators { class Advection; }
namespace Functions { class BoundaryFunction; }

namespace Energy {

class EnergyBase : public PKPhysicalBDFBase {

public:

  EnergyBase(Teuchos::Ptr<State> S, const Teuchos::RCP<Teuchos::ParameterList>& plist,
             Teuchos::ParameterList& FElist,
             const Teuchos::RCP<TreeVector>& solution);
  
  // Virtual destructor
  virtual ~EnergyBase() {}

  // EnergyBase is a PK
  // -- Setup data
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);

  // -- Calculate any diagnostics prior to doing vis
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S);


  // EnergyBase is a BDFFnBase
  // computes the non-linear functional f = f(t,u,udot)
  virtual void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f);

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // problems with temperatures -- setting a range of admissible temps
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up);

  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);
    
  // evaluating consistent faces for given BCs and cell values
  virtual void CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u);

  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du);

 protected:
  // These must be provided by the deriving PK.
  // -- setup the evaluators
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) = 0;

  // -- get enthalpy as a function of Dirichlet boundary data.  Note that this
  //    will get replaced by a better system when we get maps on the boundary
  //    faces.
  virtual void ApplyDirichletBCsToEnthalpy_(const Teuchos::Ptr<State>& S);

  // -- Add any source terms into the residual.
  virtual void AddSources_(const Teuchos::Ptr<State>& S,
                           const Teuchos::Ptr<CompositeVector>& f);
  virtual void AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h);

  // Standard methods
  virtual void SetupEnergy_(const Teuchos::Ptr<State>& S);

  // Upwinding conductivities
  virtual bool UpdateConductivityData_(const Teuchos::Ptr<State>& S);
  virtual bool UpdateConductivityDerivativeData_(const Teuchos::Ptr<State>& S);


  // boundary condition members
  virtual void UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S);

  // physical methods
  // -- accumulation of energy
  virtual void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& f);

  // -- advection of enthalpy
  virtual void AddAdvection_(const Teuchos::Ptr<State>& S,
                     const Teuchos::Ptr<CompositeVector>& f, bool negate);

  // -- diffusion of temperature
  virtual void ApplyDiffusion_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& f);

 protected:
  int niter_;

  // boundary conditions
  Teuchos::RCP<Functions::BoundaryFunction> bc_temperature_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_diff_flux_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;

  std::vector<int> bc_markers_adv_;
  std::vector<double> bc_values_adv_;
  Teuchos::RCP<Operators::BCs> bc_adv_;

  // operators
  Teuchos::RCP<Operators::Upwinding> upwinding_;
  Teuchos::RCP<Operators::Upwinding> upwinding_deriv_;

  // mathematical operators
  Teuchos::RCP<Operators::Operator> matrix_; // pc in PKPhysicalBDFBase
  Teuchos::RCP<Operators::OperatorDiffusion> matrix_diff_;
  Teuchos::RCP<Operators::OperatorAdvection> matrix_adv_;

  Teuchos::RCP<Operators::OperatorDiffusion> preconditioner_diff_;
  Teuchos::RCP<Operators::OperatorAccumulation> preconditioner_acc_;
  Teuchos::RCP<Operators::OperatorAdvection> preconditioner_adv_;
  Teuchos::RCP<Operators::Operator> lin_solver_;

  // flags and control
  bool modify_predictor_with_consistent_faces_;
  bool modify_predictor_for_freezing_;
  bool modify_correction_for_freezing_;
  bool is_source_term_;
  bool is_mass_source_term_;
  bool implicit_advection_;
  bool implicit_advection_in_pc_;
  bool precon_used_;
  bool flux_exists_;
  bool jacobian_;
  
  double T_limit_;
  
  bool coupled_to_subsurface_via_temp_;
  bool coupled_to_subsurface_via_flux_;
  bool coupled_to_surface_via_temp_;
  bool coupled_to_surface_via_flux_;

  // Keys
  Key energy_key_;
  Key enthalpy_key_;
  Key denthalpy_key_;
  Key flux_key_;
  Key energy_flux_key_;
  Key adv_energy_flux_key_;
  Key conductivity_key_;
  Key uw_conductivity_key_;
  Key dconductivity_key_;
  Key duw_conductivity_key_;
  Key de_dT_key_;
  Key source_key_;
  Key dsource_dT_key_;
  Key mass_source_key_;
  Key ss_flux_key_;

};

} // namespace Energy
} // namespace Amanzi

#endif
