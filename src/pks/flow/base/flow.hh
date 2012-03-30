/*
Flow for ATS

Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Markus Berndt (berndt@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov) 
*/

#ifndef PKS_FLOW_BASE_FLOW_HH_
#define PKS_FLOW_BASE_FLOW_HH_

#include "Teuchos_RCP.hpp"

#include "boundary-function.hh"
#include "state.hh"
#include "matrix_mfd.hh"

#include "flow_constants.hh"
#include "bdf_fn_base.hh"
#include "PK.hh"

namespace Amanzi {
namespace Flow {

class Flow : public PK, public BDFFnBase {

public:
  Flow(Teuchos::ParameterList& flow_plist, const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution);
  
  // Flow is a (virtual) PK
  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::RCP<State>& S) = 0;

  virtual double get_dT() { return dt_; }
  
  // -- transfer operators -- pointer copy
  virtual void state_to_solution(const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<TreeVector>& soln);
  virtual void solution_to_state(const Teuchos::RCP<TreeVector>& soln,
                                 const Teuchos::RCP<State>& S);

  // -- Choose a time step compatible with physics.
  virtual double get_dt() { return dt_; }

  // -- Advance from state S to state S_next at time S0.time + dt.
  virtual bool advance(double dt) = 0;

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {}

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) = 0;

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) = 0;

  // computes a norm on u-du and returns the result
  virtual double enorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) = 0;

  
protected:
  // boundary condition members
  virtual void UpdateBoundaryConditions_();
  virtual void ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature);

  // gravity members
  void AddGravityFluxesToOperator_(const Teuchos::RCP<const State>& S,
          const std::vector<WhetStone::Tensor>& K, const CompositeVector& Krel,
          const Teuchos::RCP<MatrixMFD>& matrix);
  void AddGravityFluxesToVector_(const Teuchos::RCP<const State>& S,
          const std::vector<WhetStone::Tensor>& K, const CompositeVector& Krel,
          const Teuchos::RCP<Epetra_Vector>& darcy_mass_flux);

  // control members
  void ValidateBoundaryConditions_();
  void ProcessParameterList_(const Teuchos::RCP<State>& S);

  void WriteGMVfile_(Teuchos::RCP<State> S) const;

 protected:
  int internal_tests_;  // output information
 
  Teuchos::RCP<State> S_prev_;
  Teuchos::RCP<State> S_inter_;
  Teuchos::RCP<State> S_next_;

  double dt_;
  double dt0_;

  Teuchos::ParameterList flow_plist_;

  std::vector<WhetStone::Tensor> K;  // tensor of absolute permeability

  Teuchos::RCP<MatrixMFD> matrix_;
  Teuchos::RCP<MatrixMFD> preconditioner_;

  Teuchos::RCP<BoundaryFunction> bc_pressure_;
  Teuchos::RCP<BoundaryFunction> bc_head_;
  Teuchos::RCP<BoundaryFunction> bc_flux_;
  std::vector<Operators::Matrix_bc> bc_markers_;
  std::vector<double> bc_values_;
  
  double atol_;
  double rtol_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
