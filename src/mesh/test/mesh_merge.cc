/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

#include <iostream>
#include "UnitTest++.h"
#include "Teuchos_ParameterList.hpp"

#include "AmanziComm.hh"
#include "MeshFramework.hh"
#include "MeshFrameworkFactory.hh"
#include "MeshFrameworkAudit.hh"
#include "MeshAudit.hh"
#include "MeshHelpers.hh"

// #include "Mesh.hh"
#include "Mesh_MSTK.hh"

#include "framework_meshes.hh"
#include "geometry_harnesses.hh"

using namespace Amanzi;

void
merge2D(Teuchos::RCP<Mesh>& mesh1, 
        const Teuchos::RCP<Mesh>& mesh2,
        const std::vector<int>& ids_v1,
        const std::vector<int>& ids_v2,
        const std::vector<int>& ids_e1,
        const std::vector<int>& ids_e2)
{
  auto mesh1_fw = Teuchos::rcp_dynamic_cast<AmanziMesh::Mesh_MSTK>(mesh1->getMeshFramework());
  auto mesh2_fw = Teuchos::rcp_dynamic_cast<AmanziMesh::Mesh_MSTK>(mesh2->getMeshFramework());

  auto mstk1 = mesh1_fw->getMSTK();
  auto mstk2 = mesh2_fw->getMSTK();

  int nnodes1 = MESH_Num_Vertices(mstk1);

  // add nodes of mstk2 to mstk1 
  int idv(0);
  double xyz[3];
  MVertex_ptr v1, v2(nullptr);
  MAttrib_ptr atr_v = MAttrib_New(mstk2, "vertex", POINTER, MVERTEX);

  while ((v2 = MESH_Next_Vertex(mstk2, &idv))) {
    v1 = MV_New(mstk1);

    MV_Coords(v2, xyz);
    MV_Set_Coords(v1, xyz);
    MEnt_Set_AttVal(v2, atr_v, 0, 0, (void*)v1);
  }

  // add edges of mstk2 to mstk1 
  int ival, ide(0);
  double dval;
  void *vval;
  MEdge_ptr e1, e2(nullptr);
  MAttrib_ptr atr_e = MAttrib_New(mstk2, "edge", POINTER, MEDGE);

  while ((e2 = MESH_Next_Edge(mstk2, &ide))) {
    e1 = ME_New(mstk1);

    for (int i = 0; i < 2; ++i) {
      v2 = ME_Vertex(e2, i);
      MEnt_Get_AttVal(v2, atr_v, &ival, &dval, &vval);
      ME_Set_Vertex(e1, i, (MVertex_ptr)vval);
    }
    MEnt_Set_AttVal(e2, atr_e, 0, 0, (void*)e1);
  }

  // add face of mstk2 to mstk1 
  int idf(0);
  MFace_ptr f2 = nullptr;
  MAttrib_ptr atr_f = MAttrib_New(mstk2, "face", POINTER, MFACE);

  int dirs[10];
  MFace_ptr edges[10];

  while ((f2 = MESH_Next_Face(mstk2, &idf))) {
    auto f1 = MF_New(mstk1);

    List_ptr lst_e = MF_Edges(f2, 1, 0);
    int ne = List_Num_Entries(lst_e);

    for (int i = 0; i < ne; ++i) {
      e2 = List_Entry(lst_e, i);
      MEnt_Get_AttVal(e2, atr_e, &ival, &dval, &vval);

      edges[i] = (MEdge_ptr)vval;
      dirs[i] = MF_EdgeDir(f2, e2);
    }
    MF_Set_Edges(f1, ne, edges, dirs);
    List_Delete(lst_e);
  }

  // merge nodes and edges
  for (int n = 0; n < ids_v1.size(); ++n) {
    v1 = MESH_Vertex(mstk1, ids_v1[n]);
    v2 = MESH_Vertex(mstk2, ids_v2[n]);

    MEnt_Get_AttVal(v2, atr_v, &ival, &dval, &vval);
    MVs_Merge(v1, (MVertex_ptr)vval, 0);
  }

  for (int n = 0; n < ids_e1.size(); ++n) {
    e1 = MESH_Edge(mstk1, ids_e1[n]);
    e2 = MESH_Edge(mstk2, ids_e2[n]);

    MEnt_Get_AttVal(e2, atr_e, &ival, &dval, &vval);
    MEs_Merge(e1, (MEdge_ptr)vval, 0);
  }

  mesh1_fw->rebuildAll();
  mesh1 = Teuchos::rcp(new Mesh(mesh1_fw, Teuchos::rcp(new AmanziMesh::MeshAlgorithms()), Teuchos::null));
  cacheAll(*mesh1);
}


TEST(MERGE_MESHES)
{
  auto comm = getDefaultComm();

  auto mesh1 = createStructuredUnitQuad({ AmanziMesh::Framework::MSTK }, 10, 10, comm, Teuchos::null, Teuchos::null);
  auto mesh2 = createStructuredUnitQuad({ AmanziMesh::Framework::MSTK }, 5, 10, comm, Teuchos::null, Teuchos::null);

  int nnodes1_owned = mesh1->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
  int nnodes2_owned = mesh2->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);

  int nfaces1_owned = mesh1->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  int nfaces2_owned = mesh2->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

  // deform mesh2
  for (int n = 0; n < nnodes2_owned; ++n) {
    auto xp = mesh2->getNodeCoordinate(n);
    xp[0] += 1.0;
    mesh2->setNodeCoordinate(n, xp);
  }
  mesh2->recacheGeometry();

  // collect interface nodes
  std::vector<int> ids_v1, ids_v2;
  for (int n = 0; n < nnodes1_owned; ++n) {
    if ((mesh1->getNodeCoordinate(n))[0] == 1.0) ids_v1.push_back(n);
  }
  for (int n = 0; n < nnodes2_owned; ++n) {
    if ((mesh2->getNodeCoordinate(n))[0] == 1.0) ids_v2.push_back(n);
  }

  for (int n = 0; n < ids_v1.size(); ++n) {
    auto p1 = mesh1->getNodeCoordinate(ids_v1[n]);
    auto p2 = mesh2->getNodeCoordinate(ids_v2[n]);
    AMANZI_ASSERT(norm(p1 - p2) < 1e-12);
  }

  // collect interface edges
  std::vector<int> ids_e1, ids_e2;
  for (int n = 0; n < nfaces1_owned; ++n) {
    if ((mesh1->getFaceCentroid(n))[0] == 1.0) ids_e1.push_back(n);
  }
  for (int n = 0; n < nfaces2_owned; ++n) {
    if ((mesh2->getFaceCentroid(n))[0] == 1.0) ids_e2.push_back(n);
  }

  for (int n = 0; n < ids_e1.size(); ++n) {
    auto p1 = mesh1->getFaceCentroid(ids_e1[n]);
    auto p2 = mesh2->getFaceCentroid(ids_e2[n]);
    AMANZI_ASSERT(norm(p1 - p2) < 1e-12);
  }


  merge2D(mesh1, mesh2, ids_v1, ids_v2, ids_e1, ids_e2);

  int nnodes1_new = mesh1->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
  int nfaces1_new = mesh1->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  AMANZI_ASSERT(nnodes1_owned + nnodes2_owned == nnodes1_new + 11);
  AMANZI_ASSERT(nfaces1_owned + nfaces2_owned == nfaces1_new + 10);

  testMeshAudit<MeshAudit, AmanziMesh::Mesh>(mesh1);

}
