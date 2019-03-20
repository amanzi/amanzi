/* -*-  mode: c++; indent-tabs-mode: nil -*- */
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

#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

namespace Amanzi {

typedef Teuchos::Comm<int> Comm_type;
typedef Teuchos::MpiComm<int> MpiComm_type;
typedef Teuchos::RCP<const Comm_type> Comm_ptr_type;

} // namespace Amanzi

#endif
