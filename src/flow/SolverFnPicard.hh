/*
  This is the flow component of the Amanzi code.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
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
  void Residual(const Teuchos::RCP<Vector>& u, Teuchos::RCP<Vector>& r);

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
  bool ModifyCorrection(const Teuchos::RCP<const Vector>& r, 
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
void SolverFnPicard<Vector>::Residual(const Teuchos::RCP<Vector>& u, Teuchos::RCP<Vector>& r)
{
  RPK_->AssembleMatrixMFD(*u, 0.0);
  RPK_->AssemblePreconditionerMFD(*u, 0.0, 0.0);

  RPK_->matrix()->ComputeResidual(*u, *r);
}


/* ******************************************************************
* Use linear solver. 
****************************************************************** */
template<class Vector>
void SolverFnPicard<Vector>::ApplyPreconditioner(
    const Teuchos::RCP<const Vector>& v, const Teuchos::RCP<Vector>& hv)
{
  AmanziSolvers::LinearOperatorFactory<Matrix_MFD, Epetra_Vector, Epetra_BlockMap> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix_MFD, Epetra_Vector, Epetra_BlockMap> >
      solver = factory.Create("pcg", RPK_->solvers_list, RPK_->matrix(), RPK_->preconditioner());

  solver->ApplyInverse(*v, *hv);
}


/* ******************************************************************
* Calculate relaxation factor.                                                       
****************************************************************** */
template<class Vector>
double SolverFnPicard<Vector>::ErrorNorm(
    const Teuchos::RCP<const Vector>& u, const Teuchos::RCP<const Vector>& du)
{ 
  double norm_du, norm_u;
  du->NormInf(&norm_du);
  u->NormInf(&norm_u);
  return norm_du / (1e-10 + norm_u);
}


/* ******************************************************************
* Calculate relaxation factor.                                                       
****************************************************************** */
template<class Vector>
bool SolverFnPicard<Vector>::ModifyCorrection(
    const Teuchos::RCP<const Vector>& r, 
    const Teuchos::RCP<const Vector>& u, 
    const Teuchos::RCP<Vector>& du)
{ 
  double relaxation = 1.0;

  for (int c = 0; c < u->MyLength(); c++) {
    double diff = fabs((*du)[c]);
    double uold = (*u)[c];
    double unew = uold + (*du)[c];
    double umax = std::max(fabs(unew), fabs(uold));
    if (diff > 1e-2 * umax) relaxation = std::min(relaxation, 1e-2 * umax / diff);
  }

#ifdef HAVE_MPI
  double relaxation_tmp = relaxation;
  mesh_->get_comm()->MinAll(&relaxation_tmp, &relaxation, 1);  // find the global minimum
#endif

  return relaxation;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
