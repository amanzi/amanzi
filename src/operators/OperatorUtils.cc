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
    //for (int i=0; i<bndface_inds.size(); i++) std::cout<<bndface_inds[i]<<" "; std::cout<<"\n";
    
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
    // for (int i=0; i<bndface_inds.size(); i++) std::cout<<bndface_inds[i]<<" "; std::cout<<"\n";
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
  ASSERT(tv.Data() == Teuchos::null);
  int ierr(0);
  int my_dof = 0;
  for (TreeVector::const_iterator it = tv.begin(); it != tv.end(); ++it) {
    ASSERT((*it)->Data() != Teuchos::null);
    ierr |= CopyCompositeVectorToSuperVector(map, *(*it)->Data(), sv, my_dof);
    my_dof++;            
  }
  ASSERT(!ierr);
  return ierr;
}


int CopySuperVectorToTreeVector(const SuperMap& map,const Epetra_Vector& sv,
                                TreeVector& tv)
{
  ASSERT(tv.Data() == Teuchos::null);
  int ierr(0);
  int my_dof = 0;
  for (TreeVector::iterator it = tv.begin(); it != tv.end(); ++it) {
    ASSERT((*it)->Data() != Teuchos::null);
    ierr |= CopySuperVectorToCompositeVector(map, sv, *(*it)->Data(), my_dof);
    my_dof++;            
  }
  ASSERT(!ierr);
  return ierr;
}


/* ******************************************************************
* Add super vector to tree vector, subvector-by-subvector.
****************************************************************** */
int AddSuperVectorToTreeVector(const SuperMap& map,const Epetra_Vector& sv,
                               TreeVector& tv)
{
  ASSERT(tv.Data() == Teuchos::null);
  int ierr(0);
  int my_dof = 0;
  for (TreeVector::iterator it = tv.begin();
       it != tv.end(); ++it) {
    ASSERT((*it)->Data() != Teuchos::null);
    ierr |= AddSuperVectorToCompositeVector(map, sv, *(*it)->Data(), my_dof);
    my_dof++;            
  }
  ASSERT(!ierr);
  return ierr;
}


/* ******************************************************************
* Create super map: compatibility version
****************************************************************** */
Teuchos::RCP<SuperMap> CreateSuperMap(const CompositeVectorSpace& cvs, int schema, int n_dofs)
{
  std::vector<std::string> compnames;
  std::vector<int> dofnums;
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  std::vector<Teuchos::RCP<const Epetra_Map> > ghost_maps;

  if (schema & OPERATOR_SCHEMA_DOFS_FACE) {
    ASSERT(cvs.HasComponent("face"));
    compnames.push_back("face");
    dofnums.push_back(n_dofs);
    std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> > meshmaps =
        getMaps(*cvs.Mesh(), AmanziMesh::FACE);
    maps.push_back(meshmaps.first);
    ghost_maps.push_back(meshmaps.second);
  }

  if (schema & OPERATOR_SCHEMA_DOFS_CELL) {
    ASSERT(cvs.HasComponent("cell"));
    compnames.push_back("cell");
    dofnums.push_back(n_dofs);
    std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> > meshmaps =
        getMaps(*cvs.Mesh(), AmanziMesh::CELL);
    maps.push_back(meshmaps.first);
    ghost_maps.push_back(meshmaps.second);
  }

  if (schema & OPERATOR_SCHEMA_DOFS_EDGE) {
    ASSERT(cvs.HasComponent("edge"));
    compnames.push_back("edge");
    dofnums.push_back(n_dofs);
    std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> > meshmaps =
        getMaps(*cvs.Mesh(), AmanziMesh::EDGE);
    maps.push_back(meshmaps.first);
    ghost_maps.push_back(meshmaps.second);
  }

  if (schema & OPERATOR_SCHEMA_DOFS_NODE) {
    ASSERT(cvs.HasComponent("node"));
    compnames.push_back("node");
    dofnums.push_back(n_dofs);
    std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> > meshmaps =
        getMaps(*cvs.Mesh(), AmanziMesh::NODE);
    maps.push_back(meshmaps.first);
    ghost_maps.push_back(meshmaps.second);
  }

  if (schema & OPERATOR_SCHEMA_DOFS_BNDFACE) {
    ASSERT(cvs.HasComponent("boundary_face"));
    compnames.push_back("boundary_face");
    dofnums.push_back(n_dofs);

    std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> > meshmaps =
        getMaps(*cvs.Mesh(), AmanziMesh::BOUNDARY_FACE);
    
    std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> > facemaps =
        getMaps(*cvs.Mesh(), AmanziMesh::FACE);
    // std::cout<<*meshmaps.first<<"\n";
    // std::cout<<*meshmaps.second<<"\n";
    
    // maps.push_back(meshmaps.first);
    // ghost_maps.push_back(meshmaps.second);
    // if (meshmaps.second->Comm().MyPID() == 1){
    //   std::cout<<meshmaps.first->NumMyElements()<<" "<<meshmaps.second->NumMyElements()<<"\n";
    //   for (int i=0; i<meshmaps.second->NumMyElements();i++){
    //     int f = facemaps.second->LID(meshmaps.second->GID(i));
    //     //std::cout<<"face "<<i<<" "<<f<<" "<<cvs.Mesh()->face_centroid(f)<<"\n";
    //     std::cout<<"face "<<i<<" "<<f<<" : "<<cvs.Mesh()->face_centroid(f)<<"\n";
    //   }
    // }
    //exit(0);
    // int num_boundary_faces = meshmaps.first -> NumMyElements();    
    // Teuchos::RCP<Epetra_Map>  boundary_map =  Teuchos::rcp(new Epetra_Map(-1, num_boundary_faces, 0, meshmaps.first->Comm()));

    std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> > new_bnd_map =
      CreateBoundaryMaps(cvs.Mesh(), meshmaps, facemaps);
    
    maps.push_back(new_bnd_map.first);
    ghost_maps.push_back(new_bnd_map.second);


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
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  std::vector<Teuchos::RCP<const Epetra_Map> > ghost_maps;

  for (auto it = schema.begin(); it != schema.end(); ++it) {
    compnames.push_back(schema.KindToString(it->kind));
    dofnums.push_back(it->num);
    std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> > meshmaps =
        getMaps(*cvs.Mesh(), it->kind);
    maps.push_back(meshmaps.first);
    ghost_maps.push_back(meshmaps.second);
  }

  return Teuchos::rcp(new SuperMap(cvs.Comm(), compnames, dofnums, maps, ghost_maps));
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

    for (int c = 0; c < mesh.num_entities(AmanziMesh::CELL, AmanziMesh::OWNED); ++c) {
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
****************************************************************** */
  
std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> >
CreateBoundaryMaps(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                   std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> >& bnd_maps,
                   std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> >& face_maps){

  int num_boundary_faces_owned = bnd_maps.first -> NumMyElements();

  ASSERT(num_boundary_faces_owned > 0);
  
  Teuchos::RCP<Epetra_Map>  boundary_map =  Teuchos::rcp(new Epetra_Map(-1, num_boundary_faces_owned, 0, bnd_maps.first->Comm()));

  int n_ghosted = bnd_maps.second -> NumMyElements() - num_boundary_faces_owned;
  std::vector<int> gl_id(n_ghosted), pr_id(n_ghosted), lc_id(n_ghosted);

  int total_proc = mesh->get_comm()->NumProc();
  int my_pid = mesh->get_comm()->MyPID();
  std::vector<int> min_global_id(total_proc, 0), tmp(total_proc, 0);

  tmp[my_pid] = boundary_map->GID(0);

#ifdef HAVE_MPI
  const MPI_Comm& comm = mesh->get_comm()->Comm();
  MPI_Allreduce(tmp.data(), min_global_id.data(), total_proc, MPI_INT, MPI_SUM, comm);
#endif
  
  for (int n = num_boundary_faces_owned; n < bnd_maps.second -> NumMyElements(); n++){
    int f = face_maps.second->LID(bnd_maps.second->GID(n));
    gl_id[n - num_boundary_faces_owned] = face_maps.second->GID(f);
    //if (mesh->get_comm()->MyPID() == 1) std::cout << gl_id[n - num_boundary_faces_owned]<<" "<<n - num_boundary_faces_owned<<" : "<<mesh->face_centroid(f)<<"\n";
  }

  bnd_maps.first->RemoteIDList(n_ghosted, gl_id.data(), pr_id.data(), lc_id.data());


  int n_ghosted_new = num_boundary_faces_owned;
  for (int i=0; i<n_ghosted; i++){
    if (pr_id[i] >= 0){
      n_ghosted_new++;
    }
  }
  
  //std::cout<<*face_maps.second<<"\n";

  std::vector<int> global_id_ghosted(n_ghosted_new);
  for (int i=0; i<num_boundary_faces_owned; i++)  {
    global_id_ghosted[i] = boundary_map->GID(i);  
    // int f = face_maps.second -> LID(bnd_maps.first->GID(i));
    // if (my_pid==1) std::cout<<global_id_ghosted[i]<<" "<<f<<" "<<mesh->face_centroid(f)<<"\n";
  }

  int j = num_boundary_faces_owned;
  for (int i=0; i<n_ghosted; i++){
    if (pr_id[i] >= 0){
      int proc_id = pr_id[i];
      global_id_ghosted[j] = min_global_id[proc_id] + lc_id[i];
      j++;
    }
  }
  

  Teuchos::RCP<Epetra_Map>  boundary_map_ghosted =  Teuchos::rcp(new Epetra_Map(-1, n_ghosted_new, global_id_ghosted.data(), 0, bnd_maps.first->Comm()));
  

  return std::make_pair(boundary_map, boundary_map_ghosted);
}

std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> >
CreateNonuniformMaps(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                     std::pair<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map> >& uni_maps,
                     Epetra_IntVector& num_elems, Epetra_IntVector& block_start, AmanziMesh::Entity_kind kind){

  int num_owned = mesh->num_entities(kind, AmanziMesh::OWNED);
  int num_wghosted = mesh->num_entities(kind, AmanziMesh::USED);

  int num_nonuni_owned = 0;
  for (int n=0; n<num_owned; n++) {
    block_start[n] = num_nonuni_owned;
    num_nonuni_owned += num_elems[n];    
  }

  Teuchos::RCP<Epetra_Map> nonuni_map =  Teuchos::rcp(new Epetra_Map(-1, num_nonuni_owned, 0, uni_maps.first->Comm()));

  Epetra_IntVector block_start_gl_id(block_start);
  
  for (int n=0; n<num_owned; n++){
    int loc_id = block_start[n];
    block_start_gl_id[n] = nonuni_map->GID(loc_id);
  }
 
  Epetra_Import importer(*uni_maps.second, *uni_maps.first);
  
  int* data;
  block_start_gl_id.ExtractView(&data);
  Epetra_IntVector vv1(View, *uni_maps.first, data);
  block_start_gl_id.Import(vv1, importer, Insert);

  num_elems.ExtractView(&data);
  Epetra_IntVector vv2(View, *uni_maps.first, data);
  num_elems.Import(vv2, importer, Insert);
  

  int num_nonuni_wghosted = num_nonuni_owned;
  
  for (int n = num_owned; n<num_wghosted; n++){
     block_start[n] = num_nonuni_wghosted;
     num_nonuni_wghosted += num_elems[n];
  }

  std::vector<int> global_id_ghosted(num_nonuni_wghosted);
  for (int n=0; n<num_nonuni_owned; n++){
    global_id_ghosted[n] = nonuni_map->GID(n);
  }
  for (int k=num_owned; k<num_wghosted; k++){
    for (int n=0; n<num_elems[k]; n++){
      global_id_ghosted[ block_start[k] + n ] =  block_start_gl_id[k] + n;
    }
  }

  Teuchos::RCP<Epetra_Map> nonuni_map_ghosted =  Teuchos::rcp(new Epetra_Map(-1, num_nonuni_wghosted, global_id_ghosted.data(), 0,uni_maps.first->Comm()));
                                          
  return std::make_pair(nonuni_map,  nonuni_map_ghosted);

}

}  // namespace Operators
}  // namespace Amanzi
