/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Includes and a few helper functions to make it easier to work with Comms.

#pragma once

#include "Teuchos_RCP.hpp"
#include "AmanziTypes.hh"

#ifdef HAVE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#else
#  include "Teuchos_SerialComm.hpp"
#endif


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
  return Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
}

//
// Get a serial communicator, based on MPI_COMM_SELF if possible.
// -----------------------------------------------------------------------------
inline Comm_ptr_type
getCommSelf()
{
  return Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_SELF));
}

//
// Get a communicator, based on MPI_Comm
// -----------------------------------------------------------------------------
inline Comm_ptr_type
getComm(MPI_Comm comm)
{
  return Teuchos::rcp(new Teuchos::MpiComm<int>(comm));
}


// //
// // Wraps a communicator for providing MinAll and MaxAll, as needed in Amesos2
// // solvers.
// // -----------------------------------------------------------------------------
// template <typename Comm>
// class CommWrapper {
//  public:
//   CommWrapper(const Comm& comm) : comm_(comm.Clone()) {}

//   // hacks required to allow Epetra and Tpetra to coexist
//   template <typename Return_type>
//   void MinAll(Return_type* local, Return_type* reduced, size_t count)
//   {
//     comm_->MinAll(local, reduced, count);
//   }

//   // hacks required to allow Epetra and Tpetra to coexist
//   template <typename Return_type>
//   void MaxAll(Return_type* local, Return_type* reduced, size_t count)
//   {
//     comm_->MaxAll(local, reduced, count);
//   }

//  private:
//   Teuchos::RCP<Comm> comm_;
// };


// template <>
// class CommWrapper<Comm_ptr_type> {
//  public:
//   inline CommWrapper(const Comm_ptr_type& comm) : comm_(comm) {}

//   // hacks required to allow Epetra and Tpetra to coexist
//   template <typename Return_type>
//   void MinAll(Return_type* local, Return_type* reduced, size_t count)
//   {
//     comm_->MinAll(local, reduced, count);
//   }

//   // hacks required to allow Epetra and Tpetra to coexist
//   template <typename Return_type>
//   void MaxAll(Return_type* local, Return_type* reduced, size_t count)
//   {
//     comm_->MaxAll(local, reduced, count);
//   }

//  private:
//   Comm_ptr_type comm_;
// };

// template <typename Comm>
// Teuchos::RCP<CommWrapper<Comm>>
// getCommWrapper(const Comm& comm)
// {
//   return Teuchos::rcp(new CommWrapper<Comm>(comm));
// }


bool
sameComm(const Comm_type& c1, const Comm_type& c2);


} // namespace Amanzi
