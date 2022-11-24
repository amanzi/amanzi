#include "UnitTest++.h"

#include "Epetra_MpiComm.h"
#include "Epetra_CrsMatrix.h"
#include "HYPRE_parcsr_ls.h"
#include "Ifpack_Hypre.h"
#include "Teuchos_RCP.hpp"

SUITE(SOLVERS)
{
  TEST(GIDS_OUT_OF_ORDER)
  {
    // This tests our patches of Ifpack for Hypre.
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    int numProcs = comm.NumProc();
    int myRank = comm.MyPID();

    //
    // Construct the map [0 0 1 1]
    // Note that the rows will be out of order
    //
    int elementList[2];
    elementList[0] = 2 * myRank + 1;
    elementList[1] = 2 * myRank;
    Epetra_Map ncMap(2 * numProcs, 2, elementList, 0, comm);

    //
    // Construct the diagonal matrix
    //
    Epetra_CrsMatrix accMat(Copy, ncMap, 2, true);
    int col;
    double val;
    col = 2 * myRank + 1;
    val = 2 * myRank + 2;
    accMat.InsertGlobalValues(col, 1, &val, &col);
    col = 2 * myRank;
    val = 2 * myRank + 1;
    accMat.InsertGlobalValues(col, 1, &val, &col);
    accMat.FillComplete();

    //
    // Populate the parameters
    //
    Ifpack_Hypre prec(&accMat);

    const double tol = 1e-7;
    prec.SetParameter(Solver, &HYPRE_PCGSetMaxIter, 100); // max iterations
    prec.SetParameter(Solver, &HYPRE_PCGSetTol, tol);     // conv. tolerance
    prec.SetParameter(Solver, &HYPRE_PCGSetTwoNorm, 1); // use the two norm as the stopping criteria
    prec.SetParameter(Solver, &HYPRE_PCGSetPrintLevel, 2); // print solve info
    prec.SetParameter(Solver, &HYPRE_PCGSetLogging, 1);
    prec.SetParameter(Solver, PCG);
    prec.SetParameter(false);

    //
    // Create the preconditioner (which is actually a PCG solver)
    //
    CHECK_EQUAL(prec.Compute(), 0);

    //
    // Create the initial guess and RHS
    // Note that these use non-contiguous maps
    //
    int numVec = 2;
    Epetra_MultiVector X(ncMap, numVec);
    Epetra_MultiVector KnownX(ncMap, numVec);
    KnownX.Random();
    Epetra_MultiVector B(ncMap, numVec);
    accMat.Apply(KnownX, B);

    //
    // Solve the linear system
    // It may not like that the vectors are not contiguous
    //
    CHECK_EQUAL(0, prec.ApplyInverse(B, X));

    for (int i = 0; i != KnownX.MyLength(); ++i) {
      CHECK_CLOSE(KnownX[0][i], X[0][i], tol * 10 * pow(10.0, numProcs));
    }
  }
}
