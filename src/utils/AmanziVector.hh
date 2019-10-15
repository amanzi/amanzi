/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! Includes and a few helper functions and full definitions of Maps.

/*!

  Include this instead of directly including Epetra_MultiVector.h to ease
  transition between linear algebra packages.

*/


#ifndef AMANZI_VECTOR_HH_
#define AMANZI_VECTOR_HH_

#include "Teuchos_RCP.hpp"
#include "AmanziTypes.hh"

#ifdef TRILINOS_TPETRA_STACK

#  include "Tpetra_Vector.hpp"
#  include "Tpetra_MultiVector.hpp"

#else // Epetra stack

#  include "Epetra_Vector.h"
#  include "Epetra_MultiVector.h"
#  include "Epetra_IntVector.h"
//#include "Epetra_MultiIntVector.h" // defined in trilinos 12.??

#endif // trilinos stack

#endif
