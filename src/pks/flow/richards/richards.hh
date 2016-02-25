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
#include "upwinding.hh"

#include "OperatorDiffusionFactory.hh"
#include "OperatorAccumulation.hh"

// #include "PK_Factory.hh"
// #include "PK_PhysicalBDF.hh"
// #include "PK_PhysicalBDF_ATS.hh"
#include "pk_factory_ats.hh"
#include "pk_physical_bdf_base.hh"

namespace Amanzi {

// forward declarations
class MPCSubsurface;
class PredictorDelegateBCFlux;
namespace WhetStone { class Tensor; }

namespace Flow {

//class Richards : public PK_PhysicalBDF_ATS {
class Richards : public PKPhysicalBDFBase {

public:
  Richards(const Teuchos::RCP<Teuchos::ParameterList>& plist,
           Teuchos::ParameterList& FElist,
           const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~Richards() {}

  // main methods

  virtual void setup(const Teuchos::Ptr<State>& S){Setup(S);};
  virtual void initialize(const Teuchos::Ptr<State>& S){Initialize(S);};
  // virtual void State_to_Solution(const Teuchos::RCP<State>& S,
  //                                TreeVector& soln){state_to_solution(S, soln);};
  // virtual void Solution_to_State(TreeVector& soln,
  //                                const Teuchos::RCP<State>& S){solution_to_state(soln, S);};
  virtual bool advance(double dt){ return PKBDFBase::advance(dt);};
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {CommitStep(0, dt, S);};

  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {CalculateDiagnostics(S);};



  // -- Setup data.
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

  // -- Update diagnostics for vis.
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S);

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  virtual void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);

  // problems with pressures -- setting a range of admissible pressures
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up);

  // evaluating consistent faces for given BCs and cell values
  virtual void CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u);

protected:
  // Create of physical evaluators.
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);
  virtual void SetupRichardsFlow_(const Teuchos::Ptr<State>& S);

  // boundary condition members
  virtual void UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S, bool kr=true);

  // -- builds tensor K, along with faced-based Krel if needed by the rel-perm method
  virtual void SetAbsolutePermeabilityTensor_(const Teuchos::Ptr<State>& S);
  virtual bool UpdatePermeabilityData_(const Teuchos::Ptr<State>& S);
  virtual bool UpdatePermeabilityDerivativeData_(const Teuchos::Ptr<State>& S);

  virtual void UpdateVelocity_(const Teuchos::Ptr<State>& S);

  // physical methods
  // -- diffusion term
  virtual void ApplyDiffusion_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& g);

  // virtual void AddVaporDiffusionResidual_(const Teuchos::Ptr<State>& S,
  //         const Teuchos::Ptr<CompositeVector>& g);
  // virtual void ComputeVaporDiffusionCoef(const Teuchos::Ptr<State>& S, 
  //                                        Teuchos::RCP<CompositeVector>& vapor_diff, 
  //                                        std::string var_name);
 


  // -- accumulation term
  virtual void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g);

  // -- Add any source terms into the residual.
  virtual void AddSources_(const Teuchos::Ptr<State>& S,
                           const Teuchos::Ptr<CompositeVector>& f);
  virtual void AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h);
  
  // -- gravity contributions to matrix or vector
  // virtual void AddGravityFluxes_(const Teuchos::Ptr<const Epetra_Vector>& g_vec,
  //         const Teuchos::Ptr<const CompositeVector>& rel_perm,
  //         const Teuchos::Ptr<const CompositeVector>& rho,
  //         const Teuchos::Ptr<Operators::MatrixMFD>& matrix);

  // virtual void AddGravityFluxes_FV_(const Teuchos::Ptr<const Epetra_Vector>& g_vec,
  //         const Teuchos::Ptr<const CompositeVector>& rel_perm,
  //         const Teuchos::Ptr<const CompositeVector>& rho,
  //         const Teuchos::Ptr<Operators::Matrix_TPFA>& matrix);

  // virtual void AddGravityFluxesToVector_(const Teuchos::Ptr<const Epetra_Vector>& g_vec,
  //         const Teuchos::Ptr<const CompositeVector>& rel_perm,
  //         const Teuchos::Ptr<const CompositeVector>& rho,
  //         const Teuchos::Ptr<CompositeVector>& mass_flux);

  // Nonlinear version of CalculateConsistentFaces()
  // virtual void CalculateConsistentFacesForInfiltration_(
  //     const Teuchos::Ptr<CompositeVector>& u);
  virtual bool ModifyPredictorConsistentFaces_(double h, Teuchos::RCP<TreeVector> u);
  virtual bool ModifyPredictorWC_(double h, Teuchos::RCP<TreeVector> u);
  virtual bool ModifyPredictorFluxBCs_(double h, Teuchos::RCP<TreeVector> u);

  // virtual void PreconWC_(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // -- Possibly modify the correction before it is applied
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<TreeVector> du);




  


  
protected:
  // control switches
  Operators::UpwindMethod Krel_method_;
  bool infiltrate_only_if_unfrozen_;
  bool modify_predictor_with_consistent_faces_;
  bool modify_predictor_wc_;
  bool symmetric_;
  bool precon_wc_;
  bool is_source_term_;
  bool source_term_is_differentiable_;
  bool explicit_source_;
  bool precon_used_;
  bool clobber_surf_kr_;
  
  // coupling terms
  bool coupled_to_surface_via_head_; // surface-subsurface Dirichlet coupler
  bool coupled_to_surface_via_flux_; // surface-subsurface Neumann coupler

  // -- water coupler coupling parameters
  double surface_head_cutoff_;
  double surface_head_cutoff_alpha_;
  double surface_head_eps_;

  // permeability
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;  // absolute permeability
  Teuchos::RCP<Operators::Upwinding> upwinding_;
  Teuchos::RCP<Operators::Upwinding> upwinding_deriv_;
  Teuchos::RCP<FlowRelations::WRMPartition> wrms_;
  bool upwind_from_prev_flux_;

  // mathematical operators
  Teuchos::RCP<Operators::Operator> matrix_; // pc in PKPhysicalBDFBase
  Teuchos::RCP<Operators::OperatorDiffusionWithGravity> matrix_diff_;
  Teuchos::RCP<Operators::OperatorDiffusionWithGravity> preconditioner_diff_;
  Teuchos::RCP<Operators::OperatorDiffusionWithGravity> face_matrix_diff_;
  Teuchos::RCP<Operators::OperatorAccumulation> preconditioner_acc_;
  Teuchos::RCP<Operators::Operator> lin_solver_;

  // flag to do jacobian and therefore coef derivs
  bool jacobian_;
  

  // residual vector for vapor diffusion
  Teuchos::RCP<CompositeVector> res_vapor;
  // note PC is in PKPhysicalBDFBase

  // custom enorm tolerances
  double flux_tol_;

  // boundary condition data
  Teuchos::RCP<Functions::BoundaryFunction> bc_pressure_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_head_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_seepage_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_seepage_infilt_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_infiltration_;

  // delegates
  bool modify_predictor_bc_flux_;
  bool modify_predictor_first_bc_flux_;
  Teuchos::RCP<PredictorDelegateBCFlux> flux_predictor_;

  // is this a dynamic mesh problem
  bool dynamic_mesh_;

  // is vapor turned on
  bool vapor_diffusion_;

  // scale for perm
  double perm_scale_;

  // limiters
  double p_limit_;
  double patm_limit_;

  // keys
  Key mass_dens_key_;
  Key molar_dens_key_;
  Key perm_key_;
  Key coef_key_, dcoef_key_;
  Key uw_coef_key_, duw_coef_key_;
  Key flux_key_;
  Key flux_dir_key_;
  Key velocity_key_;
  Key source_key_;
  Key ss_flux_key_;

 private:
  // factory registration
  static RegisteredPKFactory_ATS<Richards> reg_;

  // Richards has a friend in couplers...
  friend class Amanzi::MPCSubsurface;

};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
