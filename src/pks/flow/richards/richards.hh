/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  A base two-phase, thermal Richard's equation with water vapor.

  License: BSD
  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#ifndef PK_FLOW_RICHARDS_HH_
#define PK_FLOW_RICHARDS_HH_

#include "boundary_function.hh"
#include "bdf_time_integrator.hh"
#include "matrix_mfd.hh"

#include "pk_factory.hh"
#include "pk_physical_bdf_base.hh"

namespace Amanzi {

// forward declarations
class MPCCoupledFlowEnergy;
class MPCDiagonalFlowEnergy;
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
           const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist,solution),
      PKPhysicalBDFBase(plist, solution),
      coupled_to_surface_via_source_(false),
      coupled_to_surface_via_head_(false),
      coupled_to_surface_via_flux_(false),
      coupled_to_surface_via_residual_(false),
      infiltrate_only_if_unfrozen_(false),
      modify_predictor_with_consistent_faces_(false),
      niter_(0) {
    // set a few parameters before setup
    plist_.set("solution key", "pressure");
  }

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

  // setting the solution as changed should also communicate faces
  virtual void changed_solution();

  virtual bool modify_predictor(double h, const Teuchos::RCP<TreeVector>& u);

  // evaluating consistent faces for given BCs and cell values
  virtual void CalculateConsistentFaces_(double h, const Teuchos::Ptr<TreeVector>& u);


protected:
  // Create of physical evaluators.
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);
  virtual void SetupRichardsFlow_(const Teuchos::Ptr<State>& S);

  // boundary condition members
  virtual void UpdateBoundaryConditions_();
  virtual void ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& pres);

  // -- builds tensor K, along with faced-based Krel if needed by the rel-perm method
  virtual void SetAbsolutePermeabilityTensor_(const Teuchos::Ptr<State>& S);
  virtual bool UpdatePermeabilityData_(const Teuchos::Ptr<State>& S);

  // physical methods
  // -- flux calculation
  virtual void UpdateFlux_(const Teuchos::RCP<State>& S);

  // -- diffusion term
  virtual void ApplyDiffusion_(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<CompositeVector>& g);

  // -- accumulation term
  virtual void AddAccumulation_(const Teuchos::RCP<CompositeVector>& g);

  // -- gravity contributions to matrix or vector
  virtual void AddGravityFluxes_(const Teuchos::RCP<const Epetra_Vector>& g_vec,
          const Teuchos::RCP<const CompositeVector>& rel_perm,
          const Teuchos::RCP<const CompositeVector>& rho,
          const Teuchos::RCP<Operators::MatrixMFD>& matrix);

  virtual void AddGravityFluxesToVector_(const Teuchos::RCP<const Epetra_Vector>& g_vec,
          const Teuchos::RCP<const CompositeVector>& rel_perm,
          const Teuchos::RCP<const CompositeVector>& rho,
          const Teuchos::RCP<CompositeVector>& darcy_flux);


protected:
  // control switches
  int Krel_method_;
  bool assemble_preconditioner_;
  bool coupled_to_surface_via_source_;
  bool coupled_to_surface_via_head_;
  bool coupled_to_surface_via_flux_;
  bool coupled_to_surface_via_residual_;
  double surface_head_cutoff_;
  double surface_head_cutoff_alpha_;
  double surface_head_eps_;
  int niter_;
  bool infiltrate_only_if_unfrozen_;
  bool modify_predictor_with_consistent_faces_;


  // permeability
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;  // tensor of absolute permeability
  Teuchos::RCP<Operators::Upwinding> upwinding_;

  // mathematical operators
  Teuchos::RCP<Operators::MatrixMFD> matrix_;
  Teuchos::RCP<Operators::MatrixMFD> preconditioner_;
  double mass_atol_;
  double mass_rtol_;

  // boundary condition data
  Teuchos::RCP<Functions::BoundaryFunction> bc_pressure_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_head_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_seepage_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_infiltration_;
  std::vector<Operators::Matrix_bc> bc_markers_;
  std::vector<double> bc_values_;

 private:
  // factory registration
  static RegisteredPKFactory<Richards> reg_;

  // Richards has a friend in couplers...
  friend class Amanzi::MPCCoupledFlowEnergy;
  friend class Amanzi::MPCDiagonalFlowEnergy;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
