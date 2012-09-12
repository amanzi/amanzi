/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  A base two-phase, thermal Richard's equation with water vapor.

  License: BSD
  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#ifndef PK_FLOW_RICHARDS_HH_
#define PK_FLOW_RICHARDS_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "composite_vector.hh"
#include "tree_vector.hh"
#include "state.hh"
#include "matrix_mfd.hh"
#include "upwinding.hh"
#include "boundary_function.hh"

#include "PK.hh"
#include "pk_factory.hh"
#include "bdf_time_integrator.hh"

namespace Amanzi {
namespace Flow {

const int FLOW_RELATIVE_PERM_CENTERED = 1;
const int FLOW_RELATIVE_PERM_UPWIND_GRAVITY = 2;
const int FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX = 3;
const int FLOW_RELATIVE_PERM_ARITHMETIC_MEAN = 4;

class Richards : public PK {

public:
  Richards() {};
  Richards(Teuchos::ParameterList& flow_plist, const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution);

  // main methods
  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::RCP<State>& S);

  // -- Choose a time step compatible with physics.
  virtual double get_dt() { return dt_; }

  // -- transfer operators -- pointer copy
  virtual void state_to_solution(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& soln);
  virtual void solution_to_state(const Teuchos::RCP<TreeVector>& soln,
          const Teuchos::RCP<State>& S);

  // -- Advance from state S to state S_next at time S0.time + dt.
  virtual bool advance(double dt);

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

  // computes a norm on u-du and returns the result
  virtual double enorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

protected:
  // Create of physical evaluators.
  virtual void SetupPhysicalEvaluators_(const Teuchos::RCP<State>& S);
  virtual void SetupRichardsFlow_(const Teuchos::RCP<State>& S);

  // boundary condition members
  virtual void UpdateBoundaryConditions_();
  virtual void ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& pres);

  // bdf needs help
  virtual void DeriveFaceValuesFromCellValues_(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<CompositeVector>& pres);

  // -- builds tensor K, along with faced-based Krel if needed by the rel-perm method
  virtual void SetAbsolutePermeabilityTensor_(const Teuchos::RCP<State>& S);
  virtual void UpdatePermeabilityData_(const Teuchos::RCP<State>& S);

  // physical methods
  // -- diffusion term
  virtual void ApplyDiffusion_(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<CompositeVector>& g);

  // -- accumulation term
  virtual void AddAccumulation_(const Teuchos::RCP<CompositeVector>& g);

  // -- gravity contributions
  virtual void AddGravityFluxes_(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<Operators::MatrixMFD>& matrix);
  virtual void AddGravityFluxesToVector_(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<CompositeVector>& darcy_mass_flux);


protected:
  // control switches
  int Krel_method_;
  double dt_;
  double dt0_;

  // input parameter data
  Teuchos::ParameterList flow_plist_;

  // permeability
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;  // tensor of absolute permeability
  Teuchos::RCP<Operators::Upwinding> upwinding_;

  // mathematical operators
  Teuchos::RCP<Amanzi::BDFTimeIntegrator> time_stepper_;
  Teuchos::RCP<Operators::MatrixMFD> matrix_;
  Teuchos::RCP<Operators::MatrixMFD> preconditioner_;
  double atol_, mass_atol_;
  double rtol_, mass_rtol_;
  double time_step_reduction_factor_;

  // boundary condition data
  Teuchos::RCP<Functions::BoundaryFunction> bc_pressure_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_head_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;
  std::vector<Operators::Matrix_bc> bc_markers_;
  std::vector<double> bc_values_;

  // factory registration
  static RegisteredPKFactory<Richards> reg_;

  friend class MPCDiagonalFlowEnergy;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
