/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! Includes and a few helper functions to make it easier to work with Comms.

/*!

  Basic functionality that standardizes Epetra_Comm and Teuchos::Comm

*/


#ifndef AMANZI_COMM_HH_
#define AMANZI_COMM_HH_

#include "Teuchos_RCP.hpp"
#include "AmanziTypes.hh"


#ifdef TRILINOS_TPETRA_STACK

#  ifdef HAVE_MPI
#    include "Teuchos_DefaultMpiComm.hpp"
#  else
#    include "Teuchos_SerialComm.hpp"
#  endif

#else // Epetra stack

#  ifdef HAVE_MPI
#    include "Epetra_MpiComm.h"
#  else
#    include "Epetra_SerialComm.h"
#  endif

#endif // trilinos stack


//
// Comm Helper functions
// -----------------------------------------------------------------------------
namespace Amanzi {

//
// Get a default communicator, based on MPI_COMM_WORLD if possible.
// -----------------------------------------------------------------------------
inline Comm_ptr_type
getDefaultComm()
{
#ifdef TRILINOS_TPETRA_STACK
#  ifdef HAVE_MPI
  return Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
#  else
  return Teuchos::rcp(new Teuchos::SerialComm<int>());
#  endif
#else
#  ifdef HAVE_MPI
  return Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#  else
  return Teuchos::rcp(new Epetra_SerialComm());
#  endif
#endif
}

//
// Get a serial communicator, based on MPI_COMM_SELF if possible.
// -----------------------------------------------------------------------------
inline Comm_ptr_type
getCommSelf()
{
#ifdef TRILINOS_TPETRA_STACK
#  ifdef HAVE_MPI
  return Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_SELF));
#  else
  return Teuchos::rcp(new Teuchos::SerialComm<int>());
#  endif
#else
#  ifdef HAVE_MPI
  return Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_SELF));
#  else
  return Teuchos::rcp(new Epetra_SerialComm());
#  endif
#endif
}


} // namespace Amanzi


#endif
