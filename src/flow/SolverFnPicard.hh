/*
  This is the flow component of the Amanzi code.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#ifndef AMANZI_FLOW_SOLVER_FN_PICARD_
#define AMANZI_FLOW_SOLVER_FN_PICARD_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_Vector.h"

#include "LinearOperatorFactory.hh"
#include "SolverFnBase.hh"
#include "Richards_PK.hh"

namespace Amanzi {
namespace AmanziFlow {

template<class Vector>
class SolverFnPicard : public AmanziSolvers::SolverFnBase<Vector> {
 public:
  SolverFnPicard(Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::RCP<Richards_PK> RPK) :
      mesh_(mesh), RPK_(RPK) {};
  ~SolverFnPicard() {};

  // computes the non-linear functional r = F(u)
  void Residual(const Teuchos::RCP<Vector>& u, const Teuchos::RCP<Vector>& r);

  // preconditioner toolkit
  void ApplyPreconditioner(const Teuchos::RCP<const Vector>& v,
                           const Teuchos::RCP<Vector>& hv);
  void UpdatePreconditioner(const Teuchos::RCP<const Vector>& u) {
    RPK_->preconditioner()->UpdatePreconditioner();
  }

  // error analysis
  double ErrorNorm(const Teuchos::RCP<const Vector>& u, 
                   const Teuchos::RCP<const Vector>& du);

  // allow PK to modify a correction
  ModifyCorrectionResult ModifyCorrection(const Teuchos::RCP<const Vector>& r, 
          const Teuchos::RCP<const Vector>& u, 
          const Teuchos::RCP<Vector>& du);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Richards_PK> RPK_; 
};


/* ******************************************************************
* Flow Picard residual.
****************************************************************** */
template<class Vector>
void SolverFnPicard<Vector>::Residual(const Teuchos::RCP<Vector>& u, const Teuchos::RCP<Vector>& r)
{
  RPK_->AssembleMatrixMFD(*u, 0.0);
  RPK_->AssemblePreconditionerMFD(*u, 0.0, 0.0);

  RPK_->matrix()->ComputeNegativeResidual(*u, *r);
}


/* ******************************************************************
* Use linear solver. 
****************************************************************** */
template<class Vector>
void SolverFnPicard<Vector>::ApplyPreconditioner(
    const Teuchos::RCP<const Vector>& v, const Teuchos::RCP<Vector>& hv)
{
  AmanziSolvers::LinearOperatorFactory<FlowMatrix, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<FlowMatrix, CompositeVector, CompositeVectorSpace> >
      solver = factory.Create("AztecOO", RPK_->linear_operator_list_, RPK_->matrix(), RPK_->preconditioner());

  solver->ApplyInverse(*v, *hv);
}


/* ******************************************************************
* Calculate residual error.                                                       
****************************************************************** */
template<class Vector>
double SolverFnPicard<Vector>::ErrorNorm(
    const Teuchos::RCP<const Vector>& u, const Teuchos::RCP<const Vector>& r)
{ 
  double rnorm;
  r->Norm2(&rnorm);
  return rnorm;
}


/* ******************************************************************
* Calculate relaxation factor.                                                       
****************************************************************** */
template<class Vector>
ModifyCorrectionResult SolverFnPicard<Vector>::ModifyCorrection(
    const Teuchos::RCP<const Vector>& r, 
    const Teuchos::RCP<const Vector>& u, 
    const Teuchos::RCP<Vector>& du)
{ 
  // double relaxation = CalculateRelaxationFactor(u, du);
  double relaxation = 0.01;

  du->Scale(relaxation);
  return relaxation > 0. ? AmanziSolvers::CORRECTION_MODIFIED :
      AmanziSolvers::CORRECTION_NOT_MODIFIED;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
