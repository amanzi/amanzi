/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

//! A base class for assembled preconditioners.

#pragma once

#include "Inverse.hh"

class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_Map;

namespace Amanzi {
namespace AmanziSolvers {

using Preconditioner = Inverse<Epetra_CrsMatrix,
                               Epetra_CrsMatrix,
                               Epetra_Vector,
                               Epetra_Map>;

} // namespace AmanziSolvers
} // namespace Amanzi
