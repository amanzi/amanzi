/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! Includes and a few helper functions to make it easier to work with Comms.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

  Basic functionality that standardizes Epetra_Comm and Teuchos::Comm

*/


#ifndef AMANZI_COMM_HH_
#define AMANZI_COMM_HH_

#include "Teuchos_RCP.hpp"
#include "AmanziTypes.hh"


#ifdef TRILINOS_TPETRA_STACK
#ifdef HAVE_MPI
#include "Teuchos_MpiComm.hpp"
#else
#include "Teuchos_SerialComm.hpp"
#endif
#else
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#endif

namespace Amanzi {

// used by tests!
inline
Comm_ptr_type getDefaultComm() {
#ifdef TRILINOS_TPETRA_STACK
#ifdef HAVE_MPI
  return Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
#else
  return Teuchos::rcp(new Teuchos::SerialComm<int>());
#endif
#else
#ifdef HAVE_MPI
  return Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  return Teuchos::rcp(new Epetra_SerialComm());
#endif
#endif
}

inline
Comm_ptr_type getCommSelf() {
#ifdef TRILINOS_TPETRA_STACK
#ifdef HAVE_MPI
  return Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_SELF));
#else
  return Teuchos::rcp(new Teuchos::SerialComm<int>());
#endif
#else
#ifdef HAVE_MPI
  return Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_SELF));
#else
  return Teuchos::rcp(new Epetra_SerialComm());
#endif
#endif

}


template<typename Comm>
class CommWrapper {
 public:
  CommWrapper(const Comm& comm)
      : comm_(comm.Clone()) {}

  // CommWrapper(MPI_Comm comm)
  //     : comm_(Teuchos::rcp(new Comm(comm))) {}
  
  // hacks required to allow Epetra and Tpetra to coexist
  template<typename Return_type>
  void MinAll(Return_type* local, Return_type* reduced, size_t count) {
    comm_->MinAll(local, reduced, count);
  }

  // hacks required to allow Epetra and Tpetra to coexist
  template<typename Return_type>
  void MaxAll(Return_type* local, Return_type* reduced, size_t count) {
    comm_->MaxAll(local, reduced, count);
  }

 private:
  Teuchos::RCP<Comm> comm_;
};


template<>
class CommWrapper<Comm_ptr_type> {
 public:
  inline
  CommWrapper(const Comm_ptr_type& comm)
      : comm_(comm) {}

  // CommWrapper(MPI_Comm comm)
  //     : comm_(Teuchos::rcp(new Comm_type(comm))) {}
  
  // hacks required to allow Epetra and Tpetra to coexist
  template<typename Return_type>
  void MinAll(Return_type* local, Return_type* reduced, size_t count) {
    comm_->MinAll(local, reduced, count);
  }

  // hacks required to allow Epetra and Tpetra to coexist
  template<typename Return_type>
  void MaxAll(Return_type* local, Return_type* reduced, size_t count) {
    comm_->MaxAll(local, reduced, count);
  }

 private:
  Comm_ptr_type comm_;
};

template<typename Comm>
Teuchos::RCP<CommWrapper<Comm> >
getCommWrapper(const Comm& comm) {
  return Teuchos::rcp(new CommWrapper<Comm>(comm));
}


} // namespace amanzi



#endif
