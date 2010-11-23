#include "DiffusionPrecon.hpp"
#include "Teuchos_ParameterList.hpp"

DiffusionPrecon::DiffusionPrecon(Teuchos::RCP<DiffusionMatrix> &matrix, 
				 Teuchos::ParameterList &plist, const Epetra_Map &map)
    : D(matrix), map_(map)
{
  label = strdup("DiffusionPrecon");

  Teuchos::ParameterList ml_plist = plist.sublist("ML Parameters");

  //MLprec = new ML_Epetra::MultiLevelPreconditioner(D->Sff(), ml_plist, false);
  MLprec = new ML_Epetra::MultiLevelPreconditioner(D->Sff(), false);
  //MLprec->PrintList();
}

DiffusionPrecon::~DiffusionPrecon()
{
  delete MLprec;
}

void DiffusionPrecon::Compute()
{
  MLprec->ComputePreconditioner();
}


int DiffusionPrecon::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  // WARNING: When invoked by AztecOO the arguments X and Y may be aliased:
  // possibly the same object or different views of the same underlying data.
  // To cope with this possibility we do not assign to Y until the end.

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
  //Yf = Tf;

  // BACKWARD SUBSTITUTION
  // Yc <- (Dcc)^(-1) * (Xc - Dcf * P * Yf)
  Tf = Yf;
  D->ApplyDirichletProjection(Tf);
  D->Dcf().Multiply(false, Tf, Tc);  // this should do the required parallel comm
  Tc.Update(1.0, Xc, -1.0);
  Yc.ReciprocalMultiply(1.0, D->Dcc(), Tc, 0.0);

  delete [] fvec_ptrs;

  return 0;
}

