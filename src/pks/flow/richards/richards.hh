/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  A base two-phase, thermal Richard's equation with water vapor.

  License: BSD
  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#ifndef PK_FLOW_RICHARDS_HH_
#define PK_FLOW_RICHARDS_HH_

#include "wrm_partition.hh"
#include "boundary_function.hh"
#include "matrix_mfd.hh"

#include "pk_factory.hh"
#include "pk_physical_bdf_base.hh"

namespace Amanzi {

// forward declarations
class MPCCoupledFlowEnergy;
class MPCDiagonalFlowEnergy;
class MPCSurfaceSubsurfaceDirichletCoupler;
namespace WhetStone { class Tensor; }
namespace Operators { class Upwinding; }

namespace Flow {

const int FLOW_RELATIVE_PERM_CENTERED = 1;
const int FLOW_RELATIVE_PERM_UPWIND_GRAVITY = 2;
const int FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX = 3;
const int FLOW_RELATIVE_PERM_ARITHMETIC_MEAN = 4;

class Richards : public PKPhysicalBDFBase {

public:
  Richards(Teuchos::ParameterList& plist,
           const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~Richards() {}

  // main methods
  // -- Setup data.
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S);

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

  // error monitor
  virtual double enorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du);

  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> u);

  // evaluating consistent faces for given BCs and cell values
  virtual void CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u);

  virtual void set_preconditioner(const Teuchos::RCP<Operators::Matrix> preconditioner);

protected:
  // Create of physical evaluators.
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);
  virtual void SetupRichardsFlow_(const Teuchos::Ptr<State>& S);

  // boundary condition members
  virtual void UpdateBoundaryConditions_();
  virtual void ApplyBoundaryConditions_(const Teuchos::Ptr<CompositeVector>& pres);

  // -- builds tensor K, along with faced-based Krel if needed by the rel-perm method
  virtual void SetAbsolutePermeabilityTensor_(const Teuchos::Ptr<State>& S);
  virtual bool UpdatePermeabilityData_(const Teuchos::Ptr<State>& S);

  // physical methods
  // -- diffusion term
  virtual void ApplyDiffusion_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& g);

  // -- accumulation term
  virtual void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g);

  // -- gravity contributions to matrix or vector
  virtual void AddGravityFluxes_(const Teuchos::Ptr<const Epetra_Vector>& g_vec,
          const Teuchos::Ptr<const CompositeVector>& rel_perm,
          const Teuchos::Ptr<const CompositeVector>& rho,
          const Teuchos::Ptr<Operators::MatrixMFD>& matrix);

  virtual void AddGravityFluxesToVector_(const Teuchos::Ptr<const Epetra_Vector>& g_vec,
          const Teuchos::Ptr<const CompositeVector>& rel_perm,
          const Teuchos::Ptr<const CompositeVector>& rho,
          const Teuchos::Ptr<CompositeVector>& darcy_flux);


protected:
  enum FluxUpdateMode {
    UPDATE_FLUX_ITERATION = 0,
    UPDATE_FLUX_TIMESTEP = 1,
    UPDATE_FLUX_VIS = 2,
    UPDATE_FLUX_NEVER = 3
  };

  // control switches
  FluxUpdateMode update_flux_;
  int Krel_method_;
  bool assemble_preconditioner_;
  int niter_;
  bool infiltrate_only_if_unfrozen_;
  bool modify_predictor_with_consistent_faces_;
  bool symmetric_;

  // coupling terms
  bool coupled_to_surface_via_head_; // surface-subsurface Dirichlet coupler
  bool coupled_to_surface_via_flux_; // surface-subsurface Neumann coupler
  bool coupled_to_surface_via_full_; // surface-subsurface coupler with full PC
  bool coupled_to_surface_via_residual_; // surface-subsurface water coupler,
                                         // old overland PK
  bool coupled_to_surface_via_residual_new_; // surface-subsurface water coupler,

  // -- water coupler coupling parameters
  double surface_head_cutoff_;
  double surface_head_cutoff_alpha_;
  double surface_head_eps_;

  // permeability
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;  // absolute permeability
  Teuchos::RCP<Operators::Upwinding> upwinding_;
  Teuchos::RCP<FlowRelations::WRMPartition> wrms_;

  // mathematical operators
  Teuchos::RCP<Operators::MatrixMFD> matrix_;
  Teuchos::RCP<Operators::MatrixMFD> mfd_preconditioner_;

  // custom enorm tolerances
  double mass_atol_;
  double mass_rtol_;

  // boundary condition data
  Teuchos::RCP<Functions::BoundaryFunction> bc_pressure_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_head_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_seepage_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_infiltration_;

 private:
  // factory registration
  static RegisteredPKFactory<Richards> reg_;

  // Richards has a friend in couplers...
  friend class Amanzi::MPCCoupledFlowEnergy;
  friend class Amanzi::MPCDiagonalFlowEnergy;
  friend class Amanzi::MPCSurfaceSubsurfaceDirichletCoupler;

  // is this a dynamic mesh problem
  bool dynamic_mesh_;

};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
