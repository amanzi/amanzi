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
#include "matrix_mfd.hh"
#include "pk_physical_bdf_base.hh"

namespace Amanzi {

// forward declarations
namespace Operators { class Advection; }
namespace Functions { class BoundaryFunction; }

namespace Energy {

class EnergyBase : public PKPhysicalBDFBase {

public:
  EnergyBase(Teuchos::ParameterList& plist, const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist, solution),
      PKPhysicalBDFBase(plist, solution),
      modify_predictor_with_consistent_faces_(false),
      niter_(0) {
    if (!plist_.isParameter("primary variable key"))
      plist_.set("primary variable key", "temperature");
  }

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
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}


  // EnergyBase is a BDFFnBase
  // computes the non-linear functional f = f(t,u,udot)
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

  // problems with temperatures -- setting a range of admissible temps
  virtual bool is_admissible(Teuchos::RCP<const TreeVector> up);

  // error monitor
  virtual double enorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du);

  virtual bool modify_predictor(double h, const Teuchos::RCP<TreeVector>& u);

  // evaluating consistent faces for given BCs and cell values
  virtual void CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u);

  virtual void set_preconditioner(const Teuchos::RCP<Operators::Matrix> preconditioner);

protected:
  // These must be provided by the deriving PK.
  // -- setup the evaluators
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) = 0;

  // -- get enthalpy as a function of Dirichlet boundary data.  Note that this
  //    will get replaced by a better system when we get maps on the boundary
  //    faces.
  virtual void ApplyDirichletBCsToEnthalpy_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& enth) = 0;

  // -- Add any source terms into the residual.
  virtual void AddSources_(const Teuchos::Ptr<State>& S,
                           const Teuchos::Ptr<CompositeVector>& f) = 0;
  virtual void AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h) = 0;

  // Standard methods
  virtual void SetupEnergy_(const Teuchos::Ptr<State>& S);

  // boundary condition members
  virtual void UpdateBoundaryConditions_();
  virtual void ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature);

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
  enum FluxUpdateMode {
    UPDATE_FLUX_ITERATION = 0,
    UPDATE_FLUX_TIMESTEP = 1,
    UPDATE_FLUX_VIS = 2,
    UPDATE_FLUX_NEVER = 3
  };

  int niter_;

  // boundary conditions
  Teuchos::RCP<Functions::BoundaryFunction> bc_temperature_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;

  // operators
  Teuchos::RCP<Operators::Advection> advection_;
  Teuchos::RCP<Operators::MatrixMFD> matrix_;
  Teuchos::RCP<Operators::MatrixMFD> mfd_preconditioner_;

  // flags and control
  double dT_max_;
  FluxUpdateMode update_flux_;
  bool assemble_preconditioner_;
  bool modify_predictor_with_consistent_faces_;
  bool coupled_to_surface_via_temp_;
  bool coupled_to_surface_via_flux_;

  // Keys
  Key energy_key_;
  Key cell_vol_key_;
  Key enthalpy_key_;
  Key flux_key_;
  Key energy_flux_key_;
  Key conductivity_key_;
  Key de_dT_key_;

};

} // namespace Energy
} // namespace Amanzi

#endif
