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
#include "Mesh.hh"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"

#include "CompositeVector.hh"
#include "OperatorDefs.hh"
#include "OperatorUtils.hh"
#include "Schema.hh"
#include "SuperMap.hh"
#include "TreeVector.hh"
#include "TreeVector_Utils.hh"
#include "ParallelCommunication.hh"

namespace Amanzi {
namespace Operators {


/* ******************************************************************
*                        DEPRECATED
* Copy super vector to composite vector: complex schema version.
****************************************************************** */
int
CopyCompositeVectorToSuperVector(const SuperMap& smap,
                                 const CompositeVector& cv,
                                 Epetra_Vector& sv,
                                 const Schema& schema,
                                 int block_num)
{
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    std::string name(schema.KindToString(kind));

    for (int k = 0; k < num; ++k) {
      const std::vector<int>& inds = smap.Indices(block_num, name, k);
      const Epetra_MultiVector& data = *cv.ViewComponent(name);
      for (int n = 0; n != data.MyLength(); ++n) sv[inds[n]] = data[k][n];
    }
  }

  return 0;
}


/* ******************************************************************
*                        DEPRECATED
* Copy super vector to composite vector: complex schema version
****************************************************************** */
int
CopySuperVectorToCompositeVector(const SuperMap& smap,
                                 const Epetra_Vector& sv,
                                 CompositeVector& cv,
                                 const Schema& schema,
                                 int block_num)
{
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    std::string name(schema.KindToString(kind));

    for (int k = 0; k < num; ++k) {
      const std::vector<int>& inds = smap.Indices(block_num, name, k);
      Epetra_MultiVector& data = *cv.ViewComponent(name);
      for (int n = 0; n != data.MyLength(); ++n) data[k][n] = sv[inds[n]];
    }
  }

  return 0;
}


/* ******************************************************************
* Estimate size of the matrix graph.
****************************************************************** */
unsigned int
MaxRowSize(const AmanziMesh::Mesh& mesh, int schema, unsigned int n_dofs)
{
  unsigned int row_size(0);
  int dim = mesh.getSpaceDimension();
  if (schema & OPERATOR_SCHEMA_DOFS_FACE) {
    std::size_t i = (dim == 2) ? OPERATOR_QUAD_FACES : OPERATOR_HEX_FACES;

    for (int c = 0; c < mesh.getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
         ++c) {
      i = std::max(i, mesh.getCellNumFaces(c));
    }
    row_size += 2 * i;
  }

  if (schema & OPERATOR_SCHEMA_DOFS_CELL) {
    unsigned int i = (dim == 2) ? OPERATOR_QUAD_FACES : OPERATOR_HEX_FACES;
    row_size += i + 1;
  }

  if (schema & OPERATOR_SCHEMA_DOFS_NODE) {
    unsigned int i = (dim == 2) ? OPERATOR_QUAD_NODES : OPERATOR_HEX_NODES;
    row_size += 8 * i;
  }

  if (schema & OPERATOR_SCHEMA_INDICES) { row_size += 1; }

  return row_size * n_dofs;
}


/* ******************************************************************
* Estimate size of the matrix graph: general version
****************************************************************** */
unsigned int
MaxRowSize(const AmanziMesh::Mesh& mesh, const Schema& schema)
{
  unsigned int row_size(0);
  int dim = mesh.getSpaceDimension();

  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num, ndofs(0);
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::Entity_kind::FACE) {
      ndofs = (dim == 2) ? OPERATOR_QUAD_FACES : OPERATOR_HEX_FACES;
    } else if (kind == AmanziMesh::Entity_kind::CELL) {
      ndofs = 1;
    } else if (kind == AmanziMesh::Entity_kind::NODE) {
      ndofs = (dim == 2) ? OPERATOR_QUAD_NODES : OPERATOR_HEX_NODES;
    } else if (kind == AmanziMesh::Entity_kind::EDGE) {
      ndofs = (dim == 2) ? OPERATOR_QUAD_EDGES : OPERATOR_HEX_EDGES;
    }

    row_size += ndofs * num;
  }

  return row_size;
}


/* ******************************************************************
* Generates a composite vestor space.
****************************************************************** */
Teuchos::RCP<CompositeVectorSpace>
CreateCompositeVectorSpace(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                           const std::vector<std::string>& names,
                           const std::vector<AmanziMesh::Entity_kind>& locations,
                           const std::vector<int>& num_dofs,
                           bool ghosted)
{
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(ghosted);

  std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>> mastermaps;
  std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>> ghostmaps;

  for (int i = 0; i < locations.size(); ++i) {
    Teuchos::RCP<const Epetra_BlockMap> master_mp(&mesh->getMap(locations[i], false), false);
    mastermaps[names[i]] = master_mp;
    Teuchos::RCP<const Epetra_BlockMap> ghost_mp(&mesh->getMap(locations[i], true), false);
    ghostmaps[names[i]] = ghost_mp;
  }

  cvs->SetComponents(names, locations, mastermaps, ghostmaps, num_dofs);
  return cvs;
}


Teuchos::RCP<CompositeVectorSpace>
CreateCompositeVectorSpace(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                           std::string name,
                           AmanziMesh::Entity_kind location,
                           int num_dof,
                           bool ghosted)
{
  std::vector<std::string> names(1, name);
  std::vector<AmanziMesh::Entity_kind> locations(1, location);
  std::vector<int> num_dofs(1, num_dof);

  return CreateCompositeVectorSpace(mesh, names, locations, num_dofs, ghosted);
}

} // namespace Operators
} // namespace Amanzi
