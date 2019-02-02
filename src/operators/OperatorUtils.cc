/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
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
#include "ParallelCommunication.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Convert composite vector to/from super vector.
****************************************************************** */
int CopyCompositeVectorToSuperVector(const SuperMap& smap, const CompositeVector& cv,
                                     Epetra_Vector& sv, int dofnum)
{
  if (cv.HasComponent("face") && smap.HasComponent("face")) {
    const std::vector<int>& face_inds = smap.Indices("face", dofnum);
    const Epetra_MultiVector& data = *cv.ViewComponent("face");
    for (int f = 0; f != data.MyLength(); ++f) sv[face_inds[f]] = data[0][f];
  }

  if (cv.HasComponent("cell") && smap.HasComponent("cell")) {
    const std::vector<int>& cell_inds = smap.Indices("cell", dofnum);
    const Epetra_MultiVector& data = *cv.ViewComponent("cell");
    for (int c = 0; c != data.MyLength(); ++c) sv[cell_inds[c]] = data[0][c];
  } 

  if (cv.HasComponent("edge") && smap.HasComponent("edge")) {
    const std::vector<int>& edge_inds = smap.Indices("edge", dofnum);
    const Epetra_MultiVector& data = *cv.ViewComponent("edge");
    for (int e = 0; e != data.MyLength(); ++e) sv[edge_inds[e]] = data[0][e];
  } 

  if (cv.HasComponent("node") && smap.HasComponent("node")) {
    const std::vector<int>& node_inds = smap.Indices("node", dofnum);
    const Epetra_MultiVector& data = *cv.ViewComponent("node");
    for (int v = 0; v != data.MyLength(); ++v) sv[node_inds[v]] = data[0][v];
  }

  if (cv.HasComponent("boundary_face") && smap.HasComponent("boundary_face")) {
    const std::vector<int>& bndface_inds = smap.Indices("boundary_face", dofnum);

    
    const Epetra_MultiVector& data = *cv.ViewComponent("boundary_face");
    for (int f = 0; f != data.MyLength(); ++f) sv[bndface_inds[f]] = data[0][f];
  }
  
  return 0;
}


/* ******************************************************************
* Copy super vector to composite vector, component-by-component.
****************************************************************** */
int CopySuperVectorToCompositeVector(const SuperMap& smap, const Epetra_Vector& sv,
                                     CompositeVector& cv, int dofnum)
{
  if (cv.HasComponent("face") && smap.HasComponent("face")) {
    const std::vector<int>& face_inds = smap.Indices("face", dofnum);
    Epetra_MultiVector& data = *cv.ViewComponent("face");
    for (int f = 0; f != data.MyLength(); ++f) data[0][f] = sv[face_inds[f]];
  }

  if (cv.HasComponent("cell") && smap.HasComponent("cell")) {
    const std::vector<int>& cell_inds = smap.Indices("cell", dofnum);
    Epetra_MultiVector& data = *cv.ViewComponent("cell");
    for (int c = 0; c != data.MyLength(); ++c) data[0][c] = sv[cell_inds[c]];
  } 

  if (cv.HasComponent("edge") && smap.HasComponent("edge")) {
    const std::vector<int>& edge_inds = smap.Indices("edge", dofnum);
    Epetra_MultiVector& data = *cv.ViewComponent("edge");
    for (int e = 0; e != data.MyLength(); ++e) data[0][e] = sv[edge_inds[e]];
  } 

  if (cv.HasComponent("node") && smap.HasComponent("node")) {
    const std::vector<int>& node_inds = smap.Indices("node", dofnum);
    Epetra_MultiVector& data = *cv.ViewComponent("node");
    for (int v = 0; v != data.MyLength(); ++v) data[0][v] = sv[node_inds[v]];
  }

  if (cv.HasComponent("boundary_face") && smap.HasComponent("boundary_face")) {
    const std::vector<int>& bndface_inds = smap.Indices("boundary_face", dofnum);
    Epetra_MultiVector& data = *cv.ViewComponent("boundary_face");
    for (int f = 0; f != data.MyLength(); ++f) data[0][f] = sv[bndface_inds[f]];
  } 
  return 0;
}


/* ******************************************************************
* Add super vector to composite vector, component-by-component.
****************************************************************** */
int AddSuperVectorToCompositeVector(const SuperMap& smap, const Epetra_Vector& sv,
                                    CompositeVector& cv, int dofnum)
{
  if (cv.HasComponent("face") && smap.HasComponent("face")) {
    const std::vector<int>& face_inds = smap.Indices("face", dofnum);
    Epetra_MultiVector& data = *cv.ViewComponent("face");
    for (int f = 0; f != data.MyLength(); ++f) data[0][f] += sv[face_inds[f]];
  } 

  if (cv.HasComponent("boundary_face") && smap.HasComponent("boundary_face")) {
    const std::vector<int>& bndface_inds = smap.Indices("boundary_face", dofnum);   
    const Epetra_MultiVector& data = *cv.ViewComponent("boundary_face");
    for (int f = 0; f != data.MyLength(); ++f) data[0][f] += sv[bndface_inds[f]];
  } 
  
  if (cv.HasComponent("cell") && smap.HasComponent("cell")) {
    const std::vector<int>& cell_inds = smap.Indices("cell", dofnum);
    Epetra_MultiVector& data = *cv.ViewComponent("cell");
    for (int c = 0; c != data.MyLength(); ++c) data[0][c] += sv[cell_inds[c]];
  } 

  if (cv.HasComponent("edge") && smap.HasComponent("edge")) {
    const std::vector<int>& edge_inds = smap.Indices("edge", dofnum);
    Epetra_MultiVector& data = *cv.ViewComponent("edge");
    for (int e = 0; e != data.MyLength(); ++e) data[0][e] += sv[edge_inds[e]];
  } 

  if (cv.HasComponent("node") && smap.HasComponent("node")) {
    const std::vector<int>& node_inds = smap.Indices("node", dofnum);
    Epetra_MultiVector& data = *cv.ViewComponent("node");
    for (int v = 0; v != data.MyLength(); ++v) data[0][v] += sv[node_inds[v]];
  } 
  return 0;
}


/* ******************************************************************
*                        DEPRECATED
* Copy super vector to composite vector: complex schema version.
****************************************************************** */
int CopyCompositeVectorToSuperVector(const SuperMap& smap, const CompositeVector& cv,
                                     Epetra_Vector& sv, const Schema& schema)
{
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    std::string name(schema.KindToString(it->kind));

    for (int k = 0; k < it->num; ++k) {
      const std::vector<int>& inds = smap.Indices(name, k);
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
int CopySuperVectorToCompositeVector(const SuperMap& smap, const Epetra_Vector& sv,
                                     CompositeVector& cv, const Schema& schema)
{
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    std::string name(schema.KindToString(it->kind));

    for (int k = 0; k < it->num; ++k) {
      const std::vector<int>& inds = smap.Indices(name, k);
      Epetra_MultiVector& data = *cv.ViewComponent(name);
      for (int n = 0; n != data.MyLength(); ++n) data[k][n] = sv[inds[n]];
    }
  }

  return 0;
}


/* ******************************************************************
* Nonmember: copy TreeVector to/from Super-vector
****************************************************************** */
int CopyTreeVectorToSuperVector(const SuperMap& map, const TreeVector& tv,
                                Epetra_Vector& sv)
{
  AMANZI_ASSERT(tv.Data() == Teuchos::null);
  int ierr(0);
  int my_dof = 0;
  for (TreeVector::const_iterator it = tv.begin(); it != tv.end(); ++it) {
    AMANZI_ASSERT((*it)->Data() != Teuchos::null);
    ierr |= CopyCompositeVectorToSuperVector(map, *(*it)->Data(), sv, my_dof);
    my_dof++;            
  }
  AMANZI_ASSERT(!ierr);
  return ierr;
}


int CopySuperVectorToTreeVector(const SuperMap& map,const Epetra_Vector& sv,
                                TreeVector& tv)
{
  AMANZI_ASSERT(tv.Data() == Teuchos::null);
  int ierr(0);
  int my_dof = 0;
  for (TreeVector::iterator it = tv.begin(); it != tv.end(); ++it) {
    AMANZI_ASSERT((*it)->Data() != Teuchos::null);
    ierr |= CopySuperVectorToCompositeVector(map, sv, *(*it)->Data(), my_dof);
    my_dof++;            
  }
  AMANZI_ASSERT(!ierr);
  return ierr;
}


/* ******************************************************************
* Add super vector to tree vector, subvector-by-subvector.
****************************************************************** */
int AddSuperVectorToTreeVector(const SuperMap& map,const Epetra_Vector& sv,
                               TreeVector& tv)
{
  AMANZI_ASSERT(tv.Data() == Teuchos::null);
  int ierr(0);
  int my_dof = 0;
  for (TreeVector::iterator it = tv.begin();
       it != tv.end(); ++it) {
    AMANZI_ASSERT((*it)->Data() != Teuchos::null);
    ierr |= AddSuperVectorToCompositeVector(map, sv, *(*it)->Data(), my_dof);
    my_dof++;            
  }
  AMANZI_ASSERT(!ierr);
  return ierr;
}


Teuchos::RCP<SuperMap> CreateSuperMap(const CompositeVectorSpace& cvs, int schema, int n_dofs)
{
  std::vector<std::string> compnames;
  std::vector<int> dofnums;
  std::vector<Teuchos::RCP<const Epetra_BlockMap> > maps;
  std::vector<Teuchos::RCP<const Epetra_BlockMap> > ghost_maps;

  for (auto it = cvs.begin(); it != cvs.end(); ++it) {
    compnames.push_back(*it);
    dofnums.push_back(n_dofs);
    maps.push_back(cvs.Map(*it, false));
    ghost_maps.push_back(cvs.Map(*it, true));
  }

  return Teuchos::rcp(new SuperMap(cvs.Comm(), compnames, dofnums, maps, ghost_maps));
}


/* ******************************************************************
* Create super map: general version
****************************************************************** */
Teuchos::RCP<SuperMap> CreateSuperMap(const CompositeVectorSpace& cvs, Schema& schema)
{
  std::vector<std::string> compnames;
  std::vector<int> dofnums;
  std::vector<Teuchos::RCP<const Epetra_BlockMap> > maps;
  std::vector<Teuchos::RCP<const Epetra_BlockMap> > ghost_maps;

  for (auto it = schema.begin(); it != schema.end(); ++it) {
    compnames.push_back(schema.KindToString(it->kind));
    dofnums.push_back(it->num);
    auto meshmaps = getMaps(*cvs.Mesh(), it->kind);
    maps.push_back(meshmaps.first);
    ghost_maps.push_back(meshmaps.second);
  }

  return Teuchos::rcp(new SuperMap(cvs.Comm(), compnames, dofnums, maps, ghost_maps));
}

/* ******************************************************************
* Create super map: general version
****************************************************************** */
Teuchos::RCP<SuperMap> CreateSuperMap(std::vector<Teuchos::RCP<const CompositeVectorSpace> > cvs_vec,
                                      std::vector<std::string> cvs_names)
{

  AMANZI_ASSERT(cvs_vec.size() == cvs_names.size() );
  
  std::vector<std::string> compnames;
  std::vector<int> dofnums;
  std::vector<Teuchos::RCP<const Epetra_BlockMap> > maps;
  std::vector<Teuchos::RCP<const Epetra_BlockMap> > ghost_maps;
  
  int i = 0;
  for (auto cvs : cvs_vec){
    for (auto name = cvs->begin(); name != cvs->end(); ++name){
      compnames.push_back(cvs_names[i] + "-" + *name);
      dofnums.push_back( cvs->NumVectors(*name) );
      maps.push_back( cvs->Map(*name, false) );
      ghost_maps.push_back( cvs->Map(*name, true) );
    }
    i++;
  }

  return Teuchos::rcp(new SuperMap(cvs_vec[0]->Comm(), compnames, dofnums, maps, ghost_maps));
  
}


/* ******************************************************************
* Estimate size of the matrix graph.
****************************************************************** */
unsigned int MaxRowSize(const AmanziMesh::Mesh& mesh, int schema, unsigned int n_dofs)
{
  unsigned int row_size(0);
  int dim = mesh.space_dimension();
  if (schema & OPERATOR_SCHEMA_DOFS_FACE) {
    unsigned int i = (dim == 2) ? OPERATOR_QUAD_FACES : OPERATOR_HEX_FACES;

    for (int c = 0; c < mesh.num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED); ++c) {
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

  if (schema & OPERATOR_SCHEMA_INDICES) {
    row_size += 1;
  }

  return row_size * n_dofs;
}    


/* ******************************************************************
* Estimate size of the matrix graph: general version
****************************************************************** */
unsigned int MaxRowSize(const AmanziMesh::Mesh& mesh, Schema& schema)
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
*  bnd_map - pair of master and ghost maps boundary faces (discontinuos, subset of facemaps)
*
*  Results:
*  pair of master and ghost continuous maps of boundary faces
****************************************************************** */
std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> >
CreateBoundaryMaps(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                   std::pair<Teuchos::RCP<const Epetra_BlockMap>, Teuchos::RCP<const Epetra_BlockMap> >& face_maps,
                   std::pair<Teuchos::RCP<const Epetra_BlockMap>, Teuchos::RCP<const Epetra_BlockMap> >& bnd_maps)
{
  int num_boundary_faces_owned = bnd_maps.first->NumMyElements();

  AMANZI_ASSERT(num_boundary_faces_owned > 0);
  
  Teuchos::RCP<Epetra_Map> boundary_map = Teuchos::rcp(new Epetra_Map(-1, num_boundary_faces_owned, 0, bnd_maps.first->Comm()));

  int n_ghosted = bnd_maps.second->NumMyElements() - num_boundary_faces_owned;
  std::vector<int> gl_id(n_ghosted), pr_id(n_ghosted), lc_id(n_ghosted);

  int total_proc = mesh->get_comm()->NumProc();
  int my_pid = mesh->get_comm()->MyPID();
  std::vector<int> min_global_id(total_proc, 0), tmp(total_proc, 0);

  tmp[my_pid] = boundary_map->GID(0);

#ifdef HAVE_MPI
  const MPI_Comm& comm = mesh->get_comm()->Comm();
  MPI_Allreduce(tmp.data(), min_global_id.data(), total_proc, MPI_INT, MPI_SUM, comm);
#endif
  
  for (int n = num_boundary_faces_owned; n < bnd_maps.second->NumMyElements(); n++) {
    int f = face_maps.second->LID(bnd_maps.second->GID(n));
    gl_id[n - num_boundary_faces_owned] = face_maps.second->GID(f);
  }

  bnd_maps.first->RemoteIDList(n_ghosted, gl_id.data(), pr_id.data(), lc_id.data());

  int n_ghosted_new = num_boundary_faces_owned;
  for (int i=0; i<n_ghosted; i++) {
    if (pr_id[i] >= 0) {
      n_ghosted_new++;
    }
  }
  
  std::vector<int> global_id_ghosted(n_ghosted_new);
  for (int i=0; i<num_boundary_faces_owned; i++) {
    global_id_ghosted[i] = boundary_map->GID(i);  
  }

  int j = num_boundary_faces_owned;
  for (int i=0; i<n_ghosted; i++) {
    if (pr_id[i] >= 0) {
      int proc_id = pr_id[i];
      global_id_ghosted[j] = min_global_id[proc_id] + lc_id[i];
      j++;
    }
  }
  

  Teuchos::RCP<Epetra_Map> boundary_map_ghosted = 
    Teuchos::rcp(new Epetra_Map(-1, n_ghosted_new, global_id_ghosted.data(), 0, bnd_maps.first->Comm()));

  return std::make_pair(boundary_map, boundary_map_ghosted);
}


std::pair<Teuchos::RCP<const Epetra_BlockMap>, Teuchos::RCP<const Epetra_BlockMap> >
getMaps(const AmanziMesh::Mesh& mesh, AmanziMesh::Entity_kind location) {
  switch(location) {
    case AmanziMesh::CELL:
      return std::make_pair(Teuchos::rcpFromRef(mesh.cell_map(false)),
                            Teuchos::rcpFromRef(mesh.cell_map(true)));

    case AmanziMesh::FACE:
      return std::make_pair(Teuchos::rcpFromRef(mesh.face_map(false)),
                            Teuchos::rcpFromRef(mesh.face_map(true)));

    case AmanziMesh::EDGE:
      return std::make_pair(Teuchos::rcpFromRef(mesh.edge_map(false)),
                            Teuchos::rcpFromRef(mesh.edge_map(true)));

    case AmanziMesh::NODE:
      return std::make_pair(Teuchos::rcpFromRef(mesh.node_map(false)),
                            Teuchos::rcpFromRef(mesh.node_map(true)));

    case AmanziMesh::BOUNDARY_FACE:
      return std::make_pair(Teuchos::rcpFromRef(mesh.exterior_face_map(false)),
                            Teuchos::rcpFromRef(mesh.exterior_face_map(true)));
    default:
      AMANZI_ASSERT(false);
      return std::make_pair(Teuchos::null, Teuchos::null);
  }
}

}  // namespace Operators
}  // namespace Amanzi
