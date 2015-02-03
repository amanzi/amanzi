#ifndef AMANZI_ENERGY_TEST_SOLVER_FN_BASE01_HH_
#define AMANZI_ENERGY_TEST_SOLVER_FN_BASE01_HH_

#include <math.h>

#include "CompositeVector.hh"
#include "OperatorDiffusion.hh"
#include "SolverFnBase.hh"

// PDE: -div(K(T) \nabla T) = 0
class SolverFnBase01 : public Amanzi::AmanziSolvers::SolverFnBase<Amanzi::CompositeVector> {
 public:
  SolverFnBase01(Teuchos::RCP<Amanzi::Operators::OperatorDiffusion> op,
                 Teuchos::ParameterList& prec_list) {
    op_ = op;
    prec_list_ = prec_list;
    flux = Teuchos::rcp(new Amanzi::CompositeVector(op_->DomainMap()));
  };

  void Residual(const Teuchos::RCP<Amanzi::CompositeVector>& u,
                const Teuchos::RCP<Amanzi::CompositeVector>& f) {
    // constant accumulation term
    double dT = 0.02;
    Amanzi::CompositeVector phi(*u);
    phi.PutScalar(0.2);

    op_->UpdateMatrices(flux, u);
    op_->AddAccumulationTerm(*u, phi, dT, "cell");
    op_->ApplyBCs();
    op_->ComputeNegativeResidual(*u, *f);
  }

  void ApplyPreconditioner(const Teuchos::RCP<const Amanzi::CompositeVector>& u,
                           const Teuchos::RCP<Amanzi::CompositeVector>& hu) {
    op_->ApplyInverse(*u, *hu);
  }

  double ErrorNorm(const Teuchos::RCP<const Amanzi::CompositeVector>& u,
                   const Teuchos::RCP<const Amanzi::CompositeVector>& du) {
    double norm_du, norm_u;
    du->NormInf(&norm_du);
    return norm_du;
  }

  void UpdatePreconditioner(const Teuchos::RCP<const Amanzi::CompositeVector>& up) {
    // Calculate flux using the existing matrix.
    Teuchos::RCP<Amanzi::CompositeVector> flux = Teuchos::rcp(new Amanzi::CompositeVector(*up));
    op_->UpdateFlux(*up, *flux);

    // constant accumulation term
    double dT = 0.02;
    Amanzi::CompositeVector phi(*up);
    phi.PutScalar(0.2);

    // Calculate new matrix.
    op_->UpdateMatrices(flux, up);
    op_->AddAccumulationTerm(*up, phi, dT, "cell");
    op_->ApplyBCs();

    // Assemble matrix and calculate preconditioner.
    int schema_prec_dofs = op_->schema_prec_dofs();
    op_->AssembleMatrix(schema_prec_dofs);
    op_->InitPreconditioner("Hypre AMG", prec_list_);
  }

  void ChangedSolution() {};

 private:
  Teuchos::RCP<Amanzi::Operators::OperatorDiffusion> op_;  // matrix and preconditioner
  Teuchos::ParameterList prec_list_;
  Teuchos::RCP<Amanzi::CompositeVector> flux;
};

#endif

