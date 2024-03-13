/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Includes and a few helper functions to make it easier to work with Comms.
#include "AmanziComm.hh"

namespace Amanzi {

#ifdef HAVE_MPI

#  ifdef TRILINOS_TPETRA_STACK

bool
sameComm(const Comm_type& c1, const Comm_type& c2)
{
  auto mpi_c1 = dynamic_cast<const MpiComm_type*>(&c1);
  auto mpi_c2 = dynamic_cast<const MpiComm_type*>(&c2);
  if (mpi_c1 != nullptr && mpi_c2 != nullptr) {
    int cmp;
    MPI_Comm_compare(mpi_c1->getRawMpiComm(), mpi_c2->getRawMpiComm(), &cmp);
    // here we only care that the comms include the same processes, order does
    // not matter.
    return cmp != MPI_UNEQUAL;
  }

  auto serial_c1 = dynamic_cast<const SerialComm_type*>(&c1);
  auto serial_c2 = dynamic_cast<const SerialComm_type*>(&c2);
  if (serial_c1 != nullptr && serial_c2 != nullptr) { return true; }
  return false;
}

#  else // Epetra stack

bool
sameComm(const Comm_type& c1, const Comm_type& c2)
{
  auto mpi_c1 = dynamic_cast<const MpiComm_type*>(&c1);
  auto mpi_c2 = dynamic_cast<const MpiComm_type*>(&c2);
  if (mpi_c1 != nullptr && mpi_c2 != nullptr) {
    int cmp;
    MPI_Comm_compare(mpi_c1->Comm(), mpi_c2->Comm(), &cmp);
    return cmp != MPI_UNEQUAL;
  }

  auto serial_c1 = dynamic_cast<const SerialComm_type*>(&c1);
  auto serial_c2 = dynamic_cast<const SerialComm_type*>(&c2);
  if (serial_c1 != nullptr && serial_c2 != nullptr) { return true; }
  return false;
}

#  endif
#else

bool
sameComm(const Comm_type& c1, const Comm_type& c2)
{
  return true;
}

#endif

} // namespace Amanzi
