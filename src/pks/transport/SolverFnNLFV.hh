/*
  Transport PK

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_TRANSPORT_SOLVER_FN_NLFV_HH_
#define AMANZI_TRANSPORT_SOLVER_FN_NLFV_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_Vector.h"

#include "LinearOperatorFactory.hh"
#include "SolverFnBase.hh"
#include "Transport_PK.hh"

namespace Amanzi {
namespace Transport {

template<class Vector>
class SolverFnNLFV : public AmanziSolvers::SolverFnBase<Vector> {
 public:
  SolverFnNLFV(Teuchos::RCP<const AmanziMesh::Mesh> mesh, 
               Teuchos::RCP<Transport_PK> TPK, Teuchos::RCP<Vector> b) :
      mesh_(mesh), TPK_(TPK), b_(b) 
  {
    Teuchos::RCP<const State> S = TPK_->state();
    ws = S->GetFieldData("saturation_liquid")->ViewComponent("cell", false);
    phi = S->GetFieldData("porosity")->ViewComponent("cell", false);
  }

  ~SolverFnNLFV() {};

  // computes the non-linear functional r = F(u)
  void Residual(const Teuchos::RCP<Vector>& u, const Teuchos::RCP<Vector>& r);

  // preconditioner toolkit
  int ApplyPreconditioner(const Teuchos::RCP<const Vector>& v,
                           const Teuchos::RCP<Vector>& hv);
  void UpdatePreconditioner(const Teuchos::RCP<const Vector>& u);

  // error analysis
  double ErrorNorm(const Teuchos::RCP<const Vector>& u, 
                   const Teuchos::RCP<const Vector>& du);

  // allow PK to modify a correction
  bool ModifyCorrection(const Teuchos::RCP<const Vector>& r, 
                        const Teuchos::RCP<const Vector>& u, 
                        const Teuchos::RCP<Vector>& du);

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Transport_PK> TPK_; 
  Teuchos::RCP<Vector> b_; 
  Teuchos::RCP<const Epetra_MultiVector> ws, phi;
};


/* ******************************************************************
* Nonliner residual in NLFV
****************************************************************** */
template<class Vector>
void SolverFnNLFV<Vector>::Residual(const Teuchos::RCP<Vector>& u, const Teuchos::RCP<Vector>& r)
{
  Teuchos::RCP<Dispersion> matrix = TPK_->dispersion_matrix();
  matrix->AssembleMatrix(*u);
  matrix->AddTimeDerivative(TPK_->TimeStep(), *phi, *ws);

  matrix->Apply(*u, *r);
  r->Update(-1.0, *b_, 1.0);
}


/* ******************************************************************
* Use linear solver. 
****************************************************************** */
template<class Vector>
int SolverFnNLFV<Vector>::ApplyPreconditioner(
    const Teuchos::RCP<const Vector>& v, const Teuchos::RCP<Vector>& hv)
{
  AmanziSolvers::LinearOperatorFactory<Dispersion, Epetra_Vector, Epetra_BlockMap> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Dispersion, Epetra_Vector, Epetra_BlockMap> >
     solver = factory.Create("Dispersion Solver", TPK_->solvers_list, TPK_->dispersion_matrix());

  int ierr(0);
  solver->ApplyInverse(*v, *hv);
  return ierr;
}


/* ******************************************************************
* Use linear solver. 
****************************************************************** */
template<class Vector>
void SolverFnNLFV<Vector>::UpdatePreconditioner(const Teuchos::RCP<const Vector>& u)
{
  Teuchos::RCP<Dispersion> matrix = TPK_->dispersion_matrix();
  matrix->AssembleMatrix(*u);
  matrix->AddTimeDerivative(TPK_->TimeStep(), *phi, *ws);
  matrix->UpdatePreconditioner();
}


/* ******************************************************************
* Calculate residual error.                                                       
****************************************************************** */
template<class Vector>
double SolverFnNLFV<Vector>::ErrorNorm(
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
bool SolverFnNLFV<Vector>::ModifyCorrection(
    const Teuchos::RCP<const Vector>& r, 
    const Teuchos::RCP<const Vector>& u, 
    const Teuchos::RCP<Vector>& du)
{ 
  return false;
}

}  // namespace Transport
}  // namespace Amanzi

#endif
