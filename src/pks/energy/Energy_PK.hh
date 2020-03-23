/*
  This is the energy component of the Amanzi code. 
  This is a base class for energy equations.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_ENERGY_PK_HH_
#define AMANZI_ENERGY_PK_HH_

// TPLs
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "Key.hh"
#include "Operator.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_Diffusion.hh"
#include "PK.hh"
#include "PK_BDF.hh"
#include "PK_DomainFunction.hh"
#include "PK_PhysicalBDF.hh"
#include "primary_variable_field_evaluator.hh"
#include "Tensor.hh"
#include "TreeVector.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace Energy {

class Energy_PK : public PK_PhysicalBDF {
 public:
  Energy_PK(Teuchos::ParameterList& pk_tree,
            const Teuchos::RCP<Teuchos::ParameterList>& glist,
            const Teuchos::RCP<State>& S,
            const Teuchos::RCP<TreeVector>& soln);
  virtual ~Energy_PK() {};

  // methods required by PK interface
  virtual void Setup(const Teuchos::Ptr<State>& S) override;
  virtual void Initialize(const Teuchos::Ptr<State>& S) override;
  virtual std::string name() override { return passwd_; }

  // methods required for time integration
  // -- management of the preconditioner
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> hu) override {
    return op_preconditioner_->ApplyInverse(*u->Data(), *hu->Data());
  }

  // -- check the admissibility of a solution
  //    override with the actual admissibility check
  bool IsAdmissible(Teuchos::RCP<const TreeVector> up) override {
    return true;
  }

  // -- possibly modifies the predictor that is going to be used as a
  //    starting value for the nonlinear solve in the time integrator,
  //    the time integrator will pass the predictor that is computed
  //    using extrapolation and the time step that is used to compute
  //    this predictor this function returns true if the predictor was
  //    modified, false if not
  bool ModifyPredictor(double dt, Teuchos::RCP<const TreeVector> u0,
                       Teuchos::RCP<TreeVector> u) override {
    return false;
  }

  // -- possibly modifies the correction, after the nonlinear solver (NKA)
  //    has computed it, will return true if it did change the correction,
  //    so that the nonlinear iteration can store the modified correction
  //    and pass it to NKA so that the NKA space can be updated
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double dt, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<TreeVector> du) override {
    return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }

  // -- calling this indicates that the time integration
  //    scheme is changing the value of the solution in state.
  void ChangedSolution() override {
    temperature_eval_->SetFieldAsChanged(S_.ptr());
  }

  // other methods
  bool UpdateConductivityData(const Teuchos::Ptr<State>& S);
  void UpdateSourceBoundaryData(double T0, double T1, const CompositeVector& u);
  void ComputeBCs(const CompositeVector& u);

  // access 
  virtual Teuchos::RCP<Operators::Operator>
      my_operator(const Operators::OperatorType& type) override { return op_preconditioner_; } 

  virtual Teuchos::RCP<Operators::PDE_HelperDiscretization>
      my_pde(const Operators::PDEType& type) override { return op_matrix_diff_; } 

  // -- for unit tests
  std::vector<WhetStone::Tensor>& get_K() { return K; } 
  Teuchos::RCP<PrimaryVariableFieldEvaluator>& temperature_eval() { return temperature_eval_; }

 private:
  void InitializeFields_();

 public:
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

 protected:
  int dim;

  const Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> ep_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;

  // primary field
  std::string passwd_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> temperature_eval_;

  // names of state fields 
  Key temperature_key_;
  Key energy_key_, prev_energy_key_;
  Key enthalpy_key_, conductivity_key_;
  Key darcy_flux_key_, particle_density_key_, ie_rock_key_;
  Key mol_density_liquid_key_;

  // conductivity tensor
  std::vector<WhetStone::Tensor> K; 

  // boundary conditons
  std::vector<Teuchos::RCP<PK_DomainFunction> > bc_temperature_; 
  std::vector<Teuchos::RCP<PK_DomainFunction> > bc_flux_; 
  int dirichlet_bc_faces_;

  // operators and solvers
  Teuchos::RCP<Operators::PDE_Diffusion> op_matrix_diff_, op_preconditioner_diff_;
  Teuchos::RCP<Operators::PDE_Accumulation> op_acc_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_matrix_advection_, op_preconditioner_advection_;
  Teuchos::RCP<Operators::Operator> op_matrix_, op_preconditioner_, op_advection_;
  Teuchos::RCP<Operators::BCs> op_bc_;

  bool prec_include_enthalpy_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
