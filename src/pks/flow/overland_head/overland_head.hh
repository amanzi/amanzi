/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#ifndef PK_FLOW_OVERLAND_HEAD_HH_
#define PK_FLOW_OVERLAND_HEAD_HH_

#include "boundary_function.hh"
#include "matrix_mfd.hh"

#include "pk_factory.hh"
#include "pk_physical_bdf_base.hh"

namespace Amanzi {

namespace Operators { class Upwinding; }
class MPCSurfaceSubsurfaceDirichletCoupler;


namespace Flow {

namespace FlowRelations { class OverlandConductivityModel; }

class OverlandHeadFlow : public PKPhysicalBDFBase {

public:
  OverlandHeadFlow(Teuchos::ParameterList& plist,
               const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist, solution),
      PKPhysicalBDFBase(plist, solution),
      standalone_mode_(false),
      is_source_term_(false),
      coupled_to_subsurface_via_flux_(false),
      coupled_to_subsurface_via_head_(false),
      update_flux_(UPDATE_FLUX_ITERATION) {
    plist_.set("primary variable key", "surface_pressure");
    plist_.set("domain name", "surface");
  }

  // Virtual destructor
  virtual ~OverlandHeadFlow() {}

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

  virtual void set_preconditioner(const Teuchos::RCP<Operators::Matrix> preconditioner);

  virtual void changed_solution();

  // evaluating consistent faces for given BCs and cell values
  //  virtual void CalculateConsistentFaces(double h, const Teuchos::Ptr<TreeVector>& u);

protected:
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
  void ApplyDiffusion_(const Teuchos::Ptr<State>& S,const Teuchos::Ptr<CompositeVector>& g);
  // -- accumulation term
  void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g);
  // -- source terms
  void AddSourceTerms_(const Teuchos::Ptr<CompositeVector>& g);

  void test_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

 protected:
  friend class Amanzi::MPCSurfaceSubsurfaceDirichletCoupler;

  enum FluxUpdateMode {
    UPDATE_FLUX_ITERATION = 0,
    UPDATE_FLUX_TIMESTEP = 1,
    UPDATE_FLUX_VIS = 2,
    UPDATE_FLUX_NEVER = 3
  };

  // control switches
  bool standalone_mode_; // domain mesh == surface mesh
  FluxUpdateMode update_flux_;
  bool is_source_term_;
  bool assemble_preconditioner_;
  bool modify_predictor_with_consistent_faces_;
  bool symmetric_;

  // coupling term
  bool coupled_to_subsurface_via_flux_;
  bool coupled_to_subsurface_via_head_;

  // work data space
  Teuchos::RCP<Operators::Upwinding> upwinding_;

  // mathematical operators
  Teuchos::RCP<Operators::MatrixMFD> matrix_;
  Teuchos::RCP<Operators::MatrixMFD> mfd_preconditioner_;

  // boundary condition data
  Teuchos::RCP<Functions::BoundaryFunction> bc_pressure_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_zero_gradient_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_head_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;

  // overland conductivity model
  Teuchos::RCP<FlowRelations::OverlandConductivityModel> cond_model_;

  // factory registration
  static RegisteredPKFactory<OverlandHeadFlow> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
