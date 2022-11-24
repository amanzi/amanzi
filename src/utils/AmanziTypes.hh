/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! Typedefs to make forward declarations and interfaces a bit easier.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

  Forward declarations of types for use in more generic code.

*/


#ifndef AMANZI_TYPES_HH_
#define AMANZI_TYPES_HH_

#include <memory>
#include "Teuchos_RCPDecl.hpp"


#ifdef TRILINOS_TPETRA_STACK
#  include "Teuchos_Comm.hpp"
#  include "Teuchos_MpiComm.hpp"
#  include "Teuchos_SerialComm.hpp"
#else

#  include "Epetra_MpiComm.h"
#  include "Epetra_SerialComm.h"
#endif

namespace Teuchos {
template <typename T>
inline Teuchos::RCP<T>
rcp(std::unique_ptr<T>&& in)
{
  return Teuchos::rcp(in.release());
}
} // namespace Teuchos

namespace Amanzi {

#ifdef TRILINOS_TPETRA_STACK
typedef Teuchos::Comm<int> Comm_type;
#  ifdef HAVE_MPI
typedef Teuchos::MpiComm<int> MpiComm_type;
typedef Teuchos::SerialComm<int> SerialComm_type;
#  endif
#else
typedef Epetra_Comm Comm_type;
#  ifdef HAVE_MPI
typedef Epetra_MpiComm MpiComm_type;
typedef Epetra_SerialComm SerialComm_type;
#  endif
#endif

typedef Teuchos::RCP<const Comm_type> Comm_ptr_type;

} // namespace Amanzi

#endif
