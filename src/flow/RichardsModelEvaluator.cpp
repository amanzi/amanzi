#include "RichardsModelEvaluator.hpp"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#  include "mpi.h"
#else
#  include "Epetra_SerialComm.h"
#endif // HAVE_MPI

#include "Epetra_Vector.h"



RichardsModelEvaluator::RichardsModelEvaluator(RichardsProblem *problem, 
					       Teuchos::RCP<DiffusionMatrix> &matrix,
					       Teuchos::ParameterList &plist, 
					       const Epetra_Map &map) 
  : problem_(problem),  D(matrix), map_(map)
{
  
  
  
}

void RichardsModelEvaluator::initialize(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, Teuchos::ParameterList &params)
{
  
}

// Overridden from BDF2::fnBase

void RichardsModelEvaluator::fun(const double t, const Epetra_Vector& u, 
				 const Epetra_Vector& udot, Epetra_Vector& f) 
{
  // compute F(u)
  problem_->ComputeF(u, f);
 

  Epetra_Vector *uc     = problem_->CreateCellView(u);  
  Epetra_Vector *udotc  = problem_->CreateCellView(udot);
  Epetra_Vector *fc     = problem_->CreateCellView(f);
  
  // compute S'(p)
  Epetra_Vector dS (problem_->CellMap());
  problem_->dSofP(*uc, dS);
  
  // on the cell unknowns compute f=f+dS*udotc
  fc->Multiply(1.0,dS,*udotc,1.0);
}

void RichardsModelEvaluator::precon(const Epetra_Vector& X, Epetra_Vector& Y)
{

  const int ncell = D->Mesh().count_entities(Mesh_data::CELL, OWNED);

  const Epetra_Map &cell_map = D->Mesh().cell_map(false);
  const Epetra_Map &face_map = D->Mesh().face_map(false);

  // The cell and face-based DoF are packed together into the X and Y Epetra
  // vectors: cell-based DoF in the first part, followed by the face-based DoF.
  // In addition, only the owned DoF belong to the vectors.

  // Create views Xc and Xf into the cell and face segments of X.
  double **cvec_ptrs = X.Pointers();
  double **fvec_ptrs = new double*[X.NumVectors()];
  for (int i = 0; i < X.NumVectors(); ++i) fvec_ptrs[i] = cvec_ptrs[i] + ncell;
  Epetra_MultiVector Xc(View, cell_map, cvec_ptrs, X.NumVectors());
  Epetra_MultiVector Xf(View, face_map, fvec_ptrs, X.NumVectors());

  // Create views Yc and Yf into the cell and face segments of Y.
  cvec_ptrs = Y.Pointers();
  for (int i = 0; i < X.NumVectors(); ++i) fvec_ptrs[i] = cvec_ptrs[i] + ncell;
  Epetra_MultiVector Yc(View, cell_map, cvec_ptrs, X.NumVectors());
  Epetra_MultiVector Yf(View, face_map, fvec_ptrs, X.NumVectors());

  // Temporary cell and face vectors.
  Epetra_MultiVector Tc(cell_map, X.NumVectors());
  Epetra_MultiVector Tf(face_map, X.NumVectors());

  // FORWARD ELIMINATION
  // Tf <- Xf - P (Dcf)^T (Dcc)^(-1) Xc
  Tc.ReciprocalMultiply(1.0, D->Dcc(), Xc, 0.0);
  D->Dcf().Multiply(true, Tc, Tf);  // this should do the required parallel comm
  D->ApplyDirichletProjection(Tf);
  Tf.Update(1.0, Xf, -1.0);

  // "Solve" the Schur complement system for Yf with Tf as the rhs using ML
  MLprec->ApplyInverse(Tf, Yf);

  // BACKWARD SUBSTITUTION
  // Yc <- (Dcc)^(-1) * (Xc - Dcf * P * Yf)
  Tf = Yf;
  D->ApplyDirichletProjection(Tf);
  D->Dcf().Multiply(false, Tf, Tc);  // this should do the required parallel comm
  Tc.Update(1.0, Xc, -1.0);
  Yc.ReciprocalMultiply(1.0, D->Dcc(), Tc, 0.0);

  delete [] fvec_ptrs;

}

void RichardsModelEvaluator::update_precon(const double t, const Epetra_Vector& up, const double h, int& errc)
{
  problem_->ComputePrecon(up,h);

  errc = 0;
}



double RichardsModelEvaluator::enorm(const Epetra_Vector& u, const Epetra_Vector& du)
{
  // simply use 2-norm of the difference for now
  
  Epetra_Vector udiff(u);
  
  udiff.Update(-1.0,du,1.0);

  double error;
  udiff.Norm2(&error);

  return error;
}


bool RichardsModelEvaluator::is_admissible(const Epetra_Vector& up)
{
  return true;
}
