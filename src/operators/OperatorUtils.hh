/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

*/

#ifndef AMANZI_OPERATORS_UTILS_HH_
#define AMANZI_OPERATORS_UTILS_HH_

#include "Teuchos_RCP.hpp"
#include "Schema.hh"

class Epetra_Vector;

namespace Amanzi {

class CompositeVector;
class TreeVector;

namespace Operators {

class SuperMap;
class Schema;

// -- complex schema version (deprecated?)
int
CopyCompositeVectorToSuperVector(const SuperMap& map,
                                 const CompositeVector& cv,
                                 Epetra_Vector& sv,
                                 const Schema& schema,
                                 int block_num = 0);
int
CopySuperVectorToCompositeVector(const SuperMap& map,
                                 const Epetra_Vector& sv,
                                 CompositeVector& cv,
                                 const Schema& schema,
                                 int block_num = 0);


// Estimate the max number of unknowns per row. Note this can be an
// overestimate, but shouldn't be an underestimate.
unsigned int
MaxRowSize(const AmanziMesh::Mesh& mesh, int schema, unsigned int n_dofs = 1);
unsigned int
MaxRowSize(const AmanziMesh::Mesh& mesh, const Schema& schema);

// Nonmember helper function
std::pair<Teuchos::RCP<const Epetra_BlockMap>, Teuchos::RCP<const Epetra_BlockMap>>
getMaps(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_kind location);

// Factory of composite vector spaces
Teuchos::RCP<CompositeVectorSpace>
CreateCompositeVectorSpace(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                           const std::vector<std::string>& names,
                           const std::vector<AmanziMesh::Entity_kind>& locations,
                           const std::vector<int>& num_dofs,
                           bool ghosted = false);

Teuchos::RCP<CompositeVectorSpace>
CreateCompositeVectorSpace(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                           std::string name,
                           AmanziMesh::Entity_kind location,
                           int num_dof,
                           bool ghosted = false);

} // namespace Operators
} // namespace Amanzi


#endif
