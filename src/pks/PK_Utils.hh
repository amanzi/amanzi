/*
  Process Kernels

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Miscalleneous collection of simple non-member functions.
*/

#ifndef AMANZI_PK_UTILS_HH_
#define AMANZI_PK_UTILS_HH_

#include "State.hh"

namespace Amanzi {

// Averages permeability tensor in horizontal direction.
void PKUtils_CalculatePermeabilityFactorInWell(
    const Teuchos::Ptr<State>& S, Teuchos::RCP<Epetra_Vector>& Kxy);

AmanziGeometry::Point PKUtils_EntityCoordinates(
    int id, AmanziMesh::Entity_ID kind, const AmanziMesh::Mesh& mesh);

}  // namespace Amanzi

#endif
