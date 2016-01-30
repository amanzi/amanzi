/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATORS_UTILS_HH_
#define AMANZI_OPERATORS_UTILS_HH_

#include "Teuchos_RCP.hpp"

class Epetra_Vector;

namespace Amanzi {

class CompositeVector;
class TreeVector;

namespace Operators {

class SuperMap;

// Nonmember CompositeVector to/from Super-vector
int CopyCompositeVectorToSuperVector(const SuperMap& map, const CompositeVector& cv,
                                     Epetra_Vector& sv, int dofnum = 0);

int CopySuperVectorToCompositeVector(const SuperMap& map, const Epetra_Vector& sv,
                                     CompositeVector& cv, int dofnum = 0);

int AddSuperVectorToCompositeVector(const SuperMap& map, const Epetra_Vector& sv,
                                    CompositeVector& cv, int dofnum = 0);

// Nonmember TreeVector to/from Super-vector
int CopyTreeVectorToSuperVector(const SuperMap& map, const TreeVector& cv, Epetra_Vector& sv);
int CopySuperVectorToTreeVector(const SuperMap& map, const Epetra_Vector& sv, TreeVector& cv);

int AddSuperVectorToTreeVector(const SuperMap& map, const Epetra_Vector& sv, TreeVector& cv);

// supermap factory from CV and schema
Teuchos::RCP<SuperMap> CreateSuperMap(const CompositeVectorSpace& cv, int schema, int n_dofs);

// Estimate the max number of unknowns per row. Note this can be an
// overestimate, but shouldn't be an underestimate.
unsigned int MaxRowSize(const AmanziMesh::Mesh& mesh, int schema, unsigned int n_dofs = 1);

}  // namespace Operators
}  // namespace Amanzi


#endif
