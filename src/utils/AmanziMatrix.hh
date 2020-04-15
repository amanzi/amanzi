/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Includes and a few helper functions and full definitions of Matrices

/*!

  Include this instead of directly including Tpetra_*.hpp to ease transition
  between linear algebra packages.
  
*/


#ifndef AMANZI_VECTOR_HH_
#define AMANZI_VECTOR_HH_

#include "Teuchos_RCP.hpp"
#include "AmanziTypes.hh"

#ifdef TRILINOS_TPETRA_STACK

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"

#else // Epetra stack

#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

#endif // trilinos stack

#endif
