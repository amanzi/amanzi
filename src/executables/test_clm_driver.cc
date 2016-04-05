#include <iostream>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "ats_clm_driver.hh"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,0);
  MPI_Comm mpi_comm(MPI_COMM_WORLD);

  Amanzi::ATSCLMDriver driver;
  int col_types = 0;
  int num_cols = 25;
  int num_types = 1;
  int num_rows = 15;

  driver.Initialize(mpi_comm, &col_types, num_cols, num_types);

  double T[num_rows * num_cols];
  double sl[num_rows * num_cols];
  double si[num_rows * num_cols];
  double e_flux[num_cols];
  double w_flux[num_cols];
  double zero = 0.;

  for (int i=0; i!=num_cols; ++i) {
    e_flux[i] = 0;
    w_flux[i] = 0;
    for (int j=0; j!=num_rows; ++j) {
      int index = j + i*num_rows;
      T[index] = zero;
      sl[index] = zero;
      si[index] = zero;
    }
  }

  driver.SetInitCLMData(&T[0], &sl[0], &si[0]);
  driver.SetCLMData(&e_flux[0], &w_flux[0]);
  driver.Advance(1.);
  driver.Finalize();
}


