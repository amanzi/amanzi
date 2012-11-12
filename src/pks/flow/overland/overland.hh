/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Gianmarco Manzini
         Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#ifndef PK_FLOW_OVERLAND_HH_
#define PK_FLOW_OVERLAND_HH_

#include "Teuchos_TimeMonitor.hpp"

#include "matrix_mfd.hh"
#include "boundary_function.hh"
#include "composite_vector_function.hh"
#include "bdf_time_integrator.hh"
#include "wrm.hh"

#include "pk_factory.hh"
#include "pk_physical_bdf_base.hh"

namespace Amanzi {

namespace Operators { class Upwinding; }

namespace Flow {

class OverlandFlow : public PKPhysicalBDFBase {

public:
  OverlandFlow(Teuchos::ParameterList& plist, const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist, solution),
      PKPhysicalBDFBase(plist, solution),
      standalone_mode_(false),
      is_source_term_(false),
      is_coupling_term_(false),
      coupled_to_surface_via_residual_(false),
      surface_head_eps_(0.),
      update_flux_(UPDATE_FLUX_ITERATION) {
    plist_.set("solution key", "ponded_depth");
  }

  // Virtual destructor
  virtual ~OverlandFlow() {}

  // main methods
  // -- Initialize owned (dependent) variables.
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S);

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

  // admissible update -- ensure non-negativity of ponded depth
  virtual bool is_admissible(Teuchos::RCP<const TreeVector> up);

  // modify the predictor to ensure non-negativity of ponded depth
  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> up);

  // setting the solution as changed should also communicate faces
  virtual void changed_solution();

private:
  // setup methods
  virtual void SetupOverlandFlow_(const Teuchos::Ptr<State>& S);
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

  // boundary condition members
  virtual void UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S);
  virtual void UpdateBoundaryConditionsNoElev_(const Teuchos::Ptr<State>& S);
  virtual void ApplyBoundaryConditions_(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<CompositeVector>& pres );

  // computational concerns in managing abs, rel perm
  // -- builds tensor K, along with faced-based Krel if needed by the rel-perm method
  bool UpdatePermeabilityData_(const Teuchos::Ptr<State>& S);

  // physical methods
  // -- diffusion term
  void ApplyDiffusion_(const Teuchos::RCP<State>& S,const Teuchos::RCP<CompositeVector>& g);
  // -- accumulation term
  void AddAccumulation_(const Teuchos::RCP<CompositeVector>& g);
  // -- source terms
  void AddLoadValue_(const Teuchos::RCP<CompositeVector>& g);

  // mesh creation
  void CreateMesh_(const Teuchos::Ptr<State>& S);

  void test_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

 private:
  enum FluxUpdateMode {
    UPDATE_FLUX_ITERATION = 0,
    UPDATE_FLUX_TIMESTEP = 1,
    UPDATE_FLUX_VIS = 2,
    UPDATE_FLUX_NEVER = 3
  };

  // control switches
  bool standalone_mode_; // domain mesh == surface mesh
  FluxUpdateMode update_flux_;

  // work data space
  Teuchos::RCP<Operators::Upwinding> upwinding_;
  Teuchos::RCP<Teuchos::Time> steptime_; //timer

  // Conductivity evaluator, needed for hacking BCs for upwinding.
  double manning_exp_;
  double slope_regularization_;
  bool is_source_term_;
  bool is_coupling_term_;
  bool coupled_to_surface_via_residual_;
  double surface_head_eps_;

  // mathematical operators
  Teuchos::RCP<Operators::MatrixMFD> matrix_;
  Teuchos::RCP<Operators::MatrixMFD> preconditioner_;

  // boundary condition data
  Teuchos::RCP<Functions::BoundaryFunction> bc_pressure_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_zero_gradient_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_head_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;
  std::vector<Operators::Matrix_bc> bc_markers_;
  std::vector<double> bc_values_;

  // factory registration
  static RegisteredPKFactory<OverlandFlow> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
