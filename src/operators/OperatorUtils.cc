/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

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
 * Convert composite vector to/from super vector.
 ****************************************************************** */
int
CopyCompositeVectorToSuperVector(const SuperMap& smap,
                                 const CompositeVector& cv, Epetra_Vector& sv,
                                 int block_num)
{
  for (const auto& compname : cv) {
    if (smap.HasComponent(block_num, compname)) {
      for (int dofnum = 0; dofnum != cv.getNumVectors(compname); ++dofnum) {
        const auto& inds = smap.Indices(block_num, compname, dofnum);
        const auto& data = *cv.ViewComponent(compname, false);
        for (int f = 0; f != data.getLocalLength(); ++f)
          sv[inds[f]] = data[dofnum][f];
      }
    }
  }
  return 0;
}


/* ******************************************************************
 * Copy super vector to composite vector, component-by-component.
 ****************************************************************** */
int
CopySuperVectorToCompositeVector(const SuperMap& smap, const Epetra_Vector& sv,
                                 CompositeVector& cv, int block_num)
{
  for (const auto& compname : cv) {
    if (smap.HasComponent(block_num, compname)) {
      for (int dofnum = 0; dofnum != cv.getNumVectors(compname); ++dofnum) {
        const auto& inds = smap.Indices(block_num, compname, dofnum);
        auto& data = *cv.ViewComponent(compname, false);
        for (int f = 0; f != data.getLocalLength(); ++f)
          data[dofnum][f] = sv[inds[f]];
      }
    }
  }
  return 0;
}


/* ******************************************************************
 * Add super vector to composite vector, component-by-component.
 ****************************************************************** */
int
AddSuperVectorToCompositeVector(const SuperMap& smap, const Epetra_Vector& sv,
                                CompositeVector& cv, int block_num)
{
  for (const auto& compname : cv) {
    if (smap.HasComponent(block_num, compname)) {
      for (int dofnum = 0; dofnum != cv.getNumVectors(compname); ++dofnum) {
        const auto& inds = smap.Indices(block_num, compname, dofnum);
        auto& data = *cv.ViewComponent(compname, false);
        for (int f = 0; f != data.getLocalLength(); ++f)
          data[dofnum][f] += sv[inds[f]];
      }
    }
  }
  return 0;
}


/* ******************************************************************
 *                        DEPRECATED
 * Copy super vector to composite vector: complex schema version.
 ****************************************************************** */
int
CopyCompositeVectorToSuperVector(const SuperMap& smap,
                                 const CompositeVector& cv, Epetra_Vector& sv,
                                 const Schema& schema, int block_num)
{
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    std::string name(schema.KindToString(it->kind));

    for (int k = 0; k < it->num; ++k) {
      const std::vector<int>& inds = smap.Indices(block_num, name, k);
      const Epetra_MultiVector& data = *cv.ViewComponent(name);
      for (int n = 0; n != data.getLocalLength(); ++n) sv[inds[n]] = data[k][n];
    }
  }

  return 0;
}


/* ******************************************************************
 *                        DEPRECATED
 * Copy super vector to composite vector: complex schema version
 ****************************************************************** */
int
CopySuperVectorToCompositeVector(const SuperMap& smap, const Epetra_Vector& sv,
                                 CompositeVector& cv, const Schema& schema,
                                 int block_num)
{
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    std::string name(schema.KindToString(it->kind));
    for (int k = 0; k < it->num; ++k) {
      const std::vector<int>& inds = smap.Indices(block_num, name, k);
      Epetra_MultiVector& data = *cv.ViewComponent(name);
      for (int n = 0; n != data.getLocalLength(); ++n) data[k][n] = sv[inds[n]];
    }
  }

  return 0;
}


/* ******************************************************************
 * Nonmember: copy TreeVector to/from Super-vector
 ****************************************************************** */
int
CopyTreeVectorToSuperVector(const SuperMap& map, const TreeVector& tv,
                            Epetra_Vector& sv)
{
  int ierr(0);

  if (tv.Data().get()) {
    return CopyCompositeVectorToSuperVector(map, *tv.Data(), sv, 0);
  } else {
    auto sub_tvs = collectTreeVectorLeaves_const(tv);
    int block_num = 0;
    for (const auto& sub_tv : sub_tvs) {
      ierr |=
        CopyCompositeVectorToSuperVector(map, *sub_tv->Data(), sv, block_num);
      block_num++;
    }
    return ierr;
  }
}


int
CopySuperVectorToTreeVector(const SuperMap& map, const Epetra_Vector& sv,
                            TreeVector& tv)
{
  int ierr(0);

  if (tv.Data().get()) {
    return CopySuperVectorToCompositeVector(map, sv, *tv.Data(), 0);
  } else {
    auto sub_tvs = collectTreeVectorLeaves(tv);
    int block_num = 0;
    for (const auto& sub_tv : sub_tvs) {
      ierr |=
        CopySuperVectorToCompositeVector(map, sv, *sub_tv->Data(), block_num);
      block_num++;
    }
    return ierr;
  }
}


/* ******************************************************************
 * Add super vector to tree vector, subvector-by-subvector.
 ****************************************************************** */
int
AddSuperVectorToTreeVector(const SuperMap& map, const Epetra_Vector& sv,
                           TreeVector& tv)
{
  int ierr(0);

  if (tv.Data().get()) {
    return AddSuperVectorToCompositeVector(map, sv, *tv.Data(), 0);
  } else {
    auto sub_tvs = collectTreeVectorLeaves(tv);
    int block_num = 0;
    for (const auto& sub_tv : sub_tvs) {
      ierr |=
        AddSuperVectorToCompositeVector(map, sv, *sub_tv->Data(), block_num);
      block_num++;
    }
    return ierr;
  }
}


/* ******************************************************************
 * TBW
 ****************************************************************** */
Teuchos::RCP<SuperMap>
CreateSuperMap(const CompositeVectorSpace& cvs, int schema, int n_dofs)
{
  std::vector<CompositeVectorSpace> cvss;
  for (int i = 0; i != n_dofs; ++i) cvss.push_back(cvs);
  return Teuchos::rcp(new SuperMap(cvss));
}


/* ******************************************************************
 * Create super map: general version
 ****************************************************************** */
Teuchos::RCP<SuperMap>
CreateSuperMap(const CompositeVectorSpace& cvs, Schema& schema)
{
  return createSuperMap(cvs);
}


/* ******************************************************************
 * Estimate size of the matrix graph.
 ****************************************************************** */
unsigned int
MaxRowSize(const AmanziMesh::Mesh& mesh, int schema, unsigned int n_dofs)
{
  unsigned int row_size(0);
  int dim = mesh.space_dimension();
  if (schema & OPERATOR_SCHEMA_DOFS_FACE) {
    unsigned int i = (dim == 2) ? OPERATOR_QUAD_FACES : OPERATOR_HEX_FACES;

    for (int c = 0; c < mesh.num_entities(AmanziMesh::CELL,
                                          AmanziMesh::Parallel_type::OWNED);
         ++c) {
      i = std::max(i, mesh.cell_get_num_faces(c));
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
MaxRowSize(const AmanziMesh::Mesh& mesh, Schema& schema)
{
  unsigned int row_size(0);
  int dim = mesh.space_dimension();

  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int ndofs;
    if (it->kind == AmanziMesh::FACE) {
      ndofs = (dim == 2) ? OPERATOR_QUAD_FACES : OPERATOR_HEX_FACES;
    } else if (it->kind == AmanziMesh::CELL) {
      ndofs = 1;
    } else if (it->kind == AmanziMesh::NODE) {
      ndofs = (dim == 2) ? OPERATOR_QUAD_NODES : OPERATOR_HEX_NODES;
    } else if (it->kind == AmanziMesh::EDGE) {
      ndofs = (dim == 2) ? OPERATOR_QUAD_EDGES : OPERATOR_HEX_EDGES;
    }

    row_size += ndofs * it->num;
  }

  return row_size;
}


/* ******************************************************************
 *  Create continuous boundary maps
 *
 *  Parameters:
 *  mesh - Amanzi mesh
 *  face_maps - pair of master and ghost face maps (continuous)
 *  bnd_map - pair of master and ghost maps boundary faces (discontinuos, subset
 *of facemaps)
 *
 *  Results:
 *  pair of master and ghost continuous maps of boundary faces
 ****************************************************************** */
std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map>>
CreateBoundaryMaps(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                   std::pair<Teuchos::RCP<const Epetra_BlockMap>,
                             Teuchos::RCP<const Epetra_BlockMap>>& face_maps,
                   std::pair<Teuchos::RCP<const Epetra_BlockMap>,
                             Teuchos::RCP<const Epetra_BlockMap>>& bnd_maps)
{
  int num_boundary_faces_owned = bnd_maps.first->getNodeNumElements();

  AMANZI_ASSERT(num_boundary_faces_owned > 0);

  Teuchos::RCP<Epetra_Map> boundary_map = Teuchos::rcp(
    new Epetra_Map(-1, num_boundary_faces_owned, 0, bnd_maps.first->Comm()));

  int n_ghosted =
    bnd_maps.second->getNodeNumElements() - num_boundary_faces_owned;
  std::vector<int> gl_id(n_ghosted), pr_id(n_ghosted), lc_id(n_ghosted);

  int total_proc = mesh->get_comm()->getSize();
  int my_pid = mesh->get_comm()->getRank();
  std::vector<int> min_global_id(total_proc, 0), tmp(total_proc, 0);

  tmp[my_pid] = boundary_map->getGlobalElement(0);

  Teuchos::reduceAll(*mesh->get_comm(),
                     Teuchos::REDUCE_SUM,
                     total_proc,
                     tmp.data(),
                     min_global_id.data());

  for (int n = num_boundary_faces_owned;
       n < bnd_maps.second->getNodeNumElements();
       n++) {
    int f =
      face_maps.second->getLocalElement(bnd_maps.second->getGlobalElement(n));
    gl_id[n - num_boundary_faces_owned] = face_maps.second->getGlobalElement(f);
  }

  bnd_maps.first->RemoteIDList(
    n_ghosted, gl_id.data(), pr_id.data(), lc_id.data());

  int n_ghosted_new = num_boundary_faces_owned;
  for (int i = 0; i < n_ghosted; i++) {
    if (pr_id[i] >= 0) { n_ghosted_new++; }
  }

  std::vector<int> global_id_ghosted(n_ghosted_new);
  for (int i = 0; i < num_boundary_faces_owned; i++) {
    global_id_ghosted[i] = boundary_map->getGlobalElement(i);
  }

  int j = num_boundary_faces_owned;
  for (int i = 0; i < n_ghosted; i++) {
    if (pr_id[i] >= 0) {
      int proc_id = pr_id[i];
      global_id_ghosted[j] = min_global_id[proc_id] + lc_id[i];
      j++;
    }
  }

  Teuchos::RCP<Epetra_Map> boundary_map_ghosted = Teuchos::rcp(new Epetra_Map(
    -1, n_ghosted_new, global_id_ghosted.data(), 0, bnd_maps.first->Comm()));

  return std::make_pair(boundary_map, boundary_map_ghosted);
}


/* ******************************************************************
 * Generates a composite vestor space.
 ****************************************************************** */
Teuchos::RCP<CompositeVectorSpace>
CreateCompositeVectorSpace(
  Teuchos::RCP<const AmanziMesh::Mesh> mesh,
  const std::vector<std::string>& names,
  const std::vector<AmanziMesh::Entity_kind>& locations,
  const std::vector<int>& num_dofs, bool ghosted)
{
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(ghosted);

  std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>> mastermaps;
  std::map<std::string, Teuchos::RCP<const Epetra_BlockMap>> ghostmaps;

  for (int i = 0; i < locations.size(); ++i) {
    Teuchos::RCP<const Epetra_BlockMap> master_mp(
      &mesh->map(locations[i], false), false);
    mastermaps[names[i]] = master_mp;
    Teuchos::RCP<const Epetra_BlockMap> ghost_mp(&mesh->map(locations[i], true),
                                                 false);
    ghostmaps[names[i]] = ghost_mp;
  }

  cvs->SetComponents(names, locations, mastermaps, ghostmaps, num_dofs);
  return cvs;
}


Teuchos::RCP<CompositeVectorSpace>
CreateCompositeVectorSpace(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                           std::string name, AmanziMesh::Entity_kind location,
                           int num_dof, bool ghosted)
{
  std::vector<std::string> names(1, name);
  std::vector<AmanziMesh::Entity_kind> locations(1, location);
  std::vector<int> num_dofs(1, num_dof);

  return CreateCompositeVectorSpace(mesh, names, locations, num_dofs, ghosted);
}

} // namespace Operators
} // namespace Amanzi
