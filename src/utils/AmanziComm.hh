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

namespace Amanzi {

// used by tests!
inline
Comm_ptr_type getDefaultComm() {
  return Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
}

inline
Comm_ptr_type getCommSelf() {
  return Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_SELF));
}



} // namespace amanzi



#endif
