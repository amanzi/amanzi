/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

#include <UnitTest++.h>

#include <iostream>

#include "AmanziComm.hh"

#include "Geometry.hh"
#include "MeshException.hh"
#include "MeshFramework.hh"
#include "MeshFrameworkFactory.hh"
#include "MeshFrameworkAudit.hh"

#include "framework_meshes.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

template<class MeshAudit_type, class Mesh_type>
void
testMeshAudit(const Teuchos::RCP<Mesh_type>& mesh) {
  // run MeshAudit
  MeshAudit_type audit(mesh);
  int status = audit.Verify();
  CHECK_EQUAL(0, status);
}

template<class MeshAudit_type, class Mesh_type>
void
testGeometry2x2(const Teuchos::RCP<Mesh_type>& mesh)
{
  double exp_cell_volume[4] = {0.25,0.25,0.25,0.25};
  double exp_cell_centroid[4][2] = {{0.25,0.25},
                                      {0.25,0.75},
                                      {0.75,0.25},
                                      {0.75,0.75}};
  double exp_face_area[12] = {0.5,0.5,0.5,0.5,
    0.5,0.5,0.5,0.5,
    0.5,0.5,0.5,0.5};
  double exp_face_centroid[12][2] = {{0.25,0.0},{0.5,0.25},
                                       {0.25,0.5},{0.0,0.25},
                                       {0.5,0.75},{0.25,1.0},
                                       {0.0,0.75},{0.75,0.0},
                                       {1.0,0.25},{0.75,0.5},
                                       {1.0,0.75},{0.75,1.0}};

  // run MeshAudit
  MeshAudit_type audit(mesh);
  int status = audit.Verify();
  CHECK_EQUAL(0, status);

  int ncells = mesh->getNumEntities(Entity_kind::CELL,Parallel_type::OWNED);
  int nfaces = mesh->getNumEntities(Entity_kind::FACE,Parallel_type::ALL);
  int space_dim = 2;

  for (int i = 0; i < ncells; i++) {
    auto centroid = mesh->getCellCentroid(i);

    // Search for a cell with the same centroid in the
    // expected list of centroid
    bool found = false;
    for (int j = 0; j < ncells; j++) {
      if (std::abs(exp_cell_centroid[j][0]-centroid[0]) < 1.0e-10 &&
          std::abs(exp_cell_centroid[j][1]-centroid[1]) < 1.0e-10) {
        found = true;
        CHECK_EQUAL(exp_cell_volume[j],mesh->getCellVolume(i));
        break;
      }
    }
    CHECK_EQUAL(true, found);

    Entity_ID_List cfaces;
    mesh->getCellFaces(i, cfaces);
    AmanziGeometry::Point normal_sum(0.0,0.0);
    for (int j = 0; j < cfaces.size(); j++) {
      auto normal = mesh->getFaceNormal(cfaces[j],i);
      normal_sum += normal;
    }
    double val = AmanziGeometry::norm(normal_sum);
    CHECK_CLOSE(0., val, 1.0e-20);
  }

  for (int i = 0; i < nfaces; i++) {
    AmanziGeometry::Point centroid = mesh->getFaceCentroid(i);

    bool found = false;
    for (int j = 0; j < nfaces; j++) {
      if (std::abs(exp_face_centroid[j][0]-centroid[0]) < 1.0e-10 &&
          std::abs(exp_face_centroid[j][1]-centroid[1]) < 1.0e-10) {
        found = true;
        CHECK_CLOSE(exp_face_area[j], mesh->getFaceArea(i), 1.e-20);

        // Natural normal is well-posed
        AmanziGeometry::Point natural_normal = mesh->getFaceNormal(i);

        // Check the normal with respect to each connected cell is given as the
        // natural times the orientation.
        Entity_ID_List cellids;
        mesh->getFaceCells(i,Parallel_type::ALL,cellids);

        for (int k = 0; k < cellids.size(); k++) {
          int orientation = 0;
          auto normal_wrt_cell = mesh->getFaceNormal(i, cellids[k], &orientation);
          CHECK(natural_normal * orientation == normal_wrt_cell);

          // check the cell's outward normal is indeed outward (assumes star-convex)
          AmanziGeometry::Point cellcentroid = mesh->getCellCentroid(cellids[k]);
          AmanziGeometry::Point facecentroid = mesh->getFaceCentroid(i);
          AmanziGeometry::Point outvec = facecentroid-cellcentroid;

          double dp = outvec * normal_wrt_cell;
          dp /= (AmanziGeometry::norm(outvec) * AmanziGeometry::norm(normal_wrt_cell));
          CHECK_CLOSE(1., dp, 1e-10);
        }
        break;
      }
    }
    CHECK_EQUAL(found,true);
  }
}


template<class MeshAudit_type, class Mesh_type>
void
testGeometry2x2x2(const Teuchos::RCP<Mesh_type>& mesh)
{

  double exp_cell_volume[8] = {0.125,0.125,0.125,0.125,
                               0.125,0.125,0.125,0.125};
  double exp_cell_centroid[8][3] = {{0.25,0.25,0.25},
                                    {0.75,0.25,0.25},
                                    {0.25,0.75,0.25},
                                    {0.75,0.75,0.25},
                                    {0.25,0.25,0.75},
                                    {0.75,0.25,0.75},
                                    {0.25,0.75,0.75},
                                    {0.75,0.75,0.75}};
  double exp_face_area[36] = {0.25,0.25,0.25,0.25,
                              0.25,0.25,0.25,0.25,
                              0.25,0.25,0.25,0.25,
                              0.25,0.25,0.25,0.25,
                              0.25,0.25,0.25,0.25,
                              0.25,0.25,0.25,0.25,
                              0.25,0.25,0.25,0.25,
                              0.25,0.25,0.25,0.25,
                              0.25,0.25,0.25,0.25};
  double exp_face_centroid[36][3] = {{0.0,0.25,0.25},
                                     {0.0,0.75,0.25},
                                     {0.0,0.25,0.75},
                                     {0.0,0.75,0.75},

                                     {0.5,0.25,0.25},
                                     {0.5,0.75,0.25},
                                     {0.5,0.25,0.75},
                                     {0.5,0.75,0.75},

                                     {1.0,0.25,0.25},
                                     {1.0,0.75,0.25},
                                     {1.0,0.25,0.75},
                                     {1.0,0.75,0.75},

                                     {0.25,0.0,0.25},
                                     {0.75,0.0,0.25},
                                     {0.25,0.0,0.75},
                                     {0.75,0.0,0.75},

                                     {0.25,0.5,0.25},
                                     {0.75,0.5,0.25},
                                     {0.25,0.5,0.75},
                                     {0.75,0.5,0.75},

                                     {0.25,1.0,0.25},
                                     {0.75,1.0,0.25},
                                     {0.25,1.0,0.75},
                                     {0.75,1.0,0.75},

                                     {0.25,0.25,0.0},
                                     {0.75,0.25,0.0},
                                     {0.25,0.75,0.0},
                                     {0.75,0.75,0.0},

                                     {0.25,0.25,0.5},
                                     {0.75,0.25,0.5},
                                     {0.25,0.75,0.5},
                                     {0.75,0.75,0.5},

                                     {0.25,0.25,1.0},
                                     {0.75,0.25,1.0},
                                     {0.25,0.75,1.0},
                                     {0.75,0.75,1.0},
  };


  // run MeshAudit
  MeshAudit_type audit(mesh);
  int status = audit.Verify();
  CHECK_EQUAL(0, status);

  int ncells = mesh->getNumEntities(Entity_kind::CELL,Parallel_type::OWNED);
  int nfaces = mesh->getNumEntities(Entity_kind::FACE,Parallel_type::ALL);
  int space_dim_ = 3;

  for (int i = 0; i < ncells; i++) {
    AmanziGeometry::Point centroid = mesh->getCellCentroid(i);

    // Search for a cell with the same centroid in the
    // expected list of centroid
    bool found = false;
    for (int j = 0; j < ncells; j++) {
      if (std::abs(exp_cell_centroid[j][0]-centroid[0]) < 1.0e-10 &&
          std::abs(exp_cell_centroid[j][1]-centroid[1]) < 1.0e-10 &&
          std::abs(exp_cell_centroid[j][2]-centroid[2]) < 1.0e-10) {
        found = true;
        CHECK_EQUAL(exp_cell_volume[j], mesh->getCellVolume(i));
        break;
      }
    }
    CHECK_EQUAL(found,true);

    Entity_ID_List cfaces;
    AmanziGeometry::Point normal_sum(3), normal(3);
    mesh->getCellFaces(i,cfaces);
    normal_sum.set(0.0);

    for (int j = 0; j < cfaces.size(); j++) {
      normal = mesh->getFaceNormal(cfaces[j],i);
      normal_sum += normal;
    }

    double val = L22(normal_sum);
    CHECK_CLOSE(val,0.0,1.0e-20);
  }

  for (int i = 0; i < nfaces; i++) {
    AmanziGeometry::Point centroid = mesh->getFaceCentroid(i);

    bool found = false;
    for (int j = 0; j < nfaces; j++) {
      if (std::abs(exp_face_centroid[j][0]-centroid[0]) < 1.0e-10 &&
          std::abs(exp_face_centroid[j][1]-centroid[1]) < 1.0e-10 &&
          std::abs(exp_face_centroid[j][2]-centroid[2]) < 1.0e-10) {
        found = true;
        CHECK_EQUAL(exp_face_area[j],mesh->getFaceArea(i));

        // Check the natural normal
        AmanziGeometry::Point normal = mesh->getFaceNormal(i);

        // Check the normal with respect to each connected cell
        Entity_ID_List cellids;
        mesh->getFaceCells(i,Parallel_type::ALL,cellids);

        for (int k = 0; k < cellids.size(); k++) {
          int dir;
          AmanziGeometry::Point normal_wrt_cell =
            mesh->getFaceNormal(i, cellids[k], &dir);

          AmanziGeometry::Point normal1(normal);
          normal1 *= dir;

          CHECK_ARRAY_EQUAL(&(normal1[0]),&(normal_wrt_cell[0]),space_dim_);

          AmanziGeometry::Point cellcentroid = mesh->getCellCentroid(cellids[k]);
          AmanziGeometry::Point facecentroid = mesh->getFaceCentroid(i);
          AmanziGeometry::Point outvec = facecentroid-cellcentroid;

          double dp = outvec*normal_wrt_cell;
          dp /= (norm(outvec)*norm(normal_wrt_cell));
          CHECK_CLOSE(dp,1.0,1e-10);
        }
        break;
      }
    }

    CHECK_EQUAL(found,true);
  }

  // Now deform the mesh a little and verify that the sum of the
  // outward normals of all faces of cell is still zero
  AmanziGeometry::Point ccoords = mesh->getNodeCoordinate(13); // central node

  // Lets be sure this is the central node
  CHECK_EQUAL(ccoords[0],0.5);
  CHECK_EQUAL(ccoords[1],0.5);
  CHECK_EQUAL(ccoords[2],0.5);

  // Perturb it
  ccoords.set(0.7,0.7,0.7);
  mesh->setNodeCoordinate(13,ccoords);

  // Now check the normals
  for (int i = 0; i < ncells; i++) {
    Entity_ID_List cfaces;
    AmanziGeometry::Point normal_sum(3), normal(3);

    mesh->getCellFaces(i,cfaces);
    normal_sum.set(0.0);

    for (int j = 0; j < cfaces.size(); j++) {
      normal = mesh->getFaceNormal(cfaces[j],i);
      normal_sum += normal;
    }

    double val = L22(normal_sum);
    CHECK_CLOSE(val,0.0,1.0e-20);
  }
}


template<class MeshAudit_type, class Mesh_type>
void
testGeometry1x1x1(const Teuchos::RCP<Mesh_type>& mesh)
{
  // run MeshAudit
  MeshAudit_type audit(mesh);
  int status = audit.Verify();
  CHECK_EQUAL(0, status);

  double xyz[12][3] = {{0, 0, 0},
		       {1, 0, 0},
		       {0, 1, 0},
		       {1, 1, 0},
		       {0, 0, 1},
		       {1, 0, 1},
		       {0, 1, 1},
		       {1, 1, 1}};
  Amanzi::AmanziMesh::Entity_ID local_cellnodes[8] = {0,1,2,3,4,5,6,7};
  Amanzi::AmanziMesh::Entity_ID local_facenodes[6][4] = {{0,1,5,4},
					{1,2,6,5},
					{2,3,7,6},
					{3,0,4,7},
					{0,3,2,1},
					{4,5,6,7}};

  // 1 cell meshes don't parallelize well...
  if (mesh->get_comm()->NumProc() > 1) return;
  int ncells = mesh->getNumEntities(Entity_kind::CELL,Parallel_type::OWNED);
  CHECK_EQUAL(1, ncells);
  int nfaces = mesh->getNumEntities(Entity_kind::FACE,Parallel_type::OWNED);
  CHECK_EQUAL(6, nfaces);
  int nnodes = mesh->getNumEntities(Entity_kind::NODE,Parallel_type::OWNED);
  CHECK_EQUAL(8, nnodes);

  for (int i = 0; i < nnodes; i++) {
    Amanzi::AmanziGeometry::Point coords(mesh->get_space_dimension());
    coords = mesh->getNodeCoordinate(i);
    CHECK_ARRAY_EQUAL(xyz[i],coords,3);
  }

  Entity_ID_List cellnodes;
  mesh->getCellNodes(0, cellnodes);
  auto ccoords = mesh->getCellCoordinates(0);
  for (int j = 0; j < 8; j++) {
    CHECK_ARRAY_EQUAL(xyz[cellnodes[j]],ccoords[j],3);
  }

  if (mesh->is_ordered()) {
    Entity_ID_List faces;
    Entity_Direction_List facedirs;
    mesh->getCellFacesAndDirs(0,faces,&facedirs);
    for (int j = 0; j < 6; j++) {
      Entity_ID_List facenodes;
      mesh->getFaceNodes(faces[j], facenodes);
      auto fcoords = mesh->getFaceCoordinates(faces[j]);

      Entity_ID_List expfacenodes(4);
      for (int k = 0; k < 4; k++)
        expfacenodes[k] = cellnodes[local_facenodes[j][k]];

      // The order of nodes returned may be different from what we expected
      // So make sure we have a matching node to start with
      int k0 = -1;
      int found = 0;
      for (int k = 0; k < 4; k++) {
        if (expfacenodes[k] == facenodes[0]) {
          k0 = k;
          found = 1;
          break;
        }
      }
      CHECK_EQUAL(found,1);

      if (facedirs[j] == 1) {
        for (int k = 0; k < 4; k++) {
          CHECK_EQUAL(expfacenodes[(k0+k)%4],facenodes[k]);
          CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0+k)%4]],fcoords[k],3);
        }
      } else {
        for (int k = 0; k < 4; k++) {
          CHECK_EQUAL(expfacenodes[(k0+4-k)%4],facenodes[k]);
          CHECK_ARRAY_EQUAL(xyz[expfacenodes[(k0+4-k)%4]],fcoords[k],3);
        }
      }
    }
  }
}


TEST(MESH_GEOMETRY_PLANAR)
{
  // only works in MSTK
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 2D geometry with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;

    auto mesh = createFrameworkStructuredUnitSquare(Preference{frm}, 2);
    testGeometry2x2<MeshFrameworkAudit, MeshFramework>(mesh);
  }
}


TEST(MESH_GEOMETRY_1BOX_GENERATED)
{
  // only works in MSTK
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }

  if (getDefaultComm()->NumProc() == 1)
    frameworks.push_back(Framework::SIMPLE);

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D geometry with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createFrameworkStructuredUnitCube(Preference{frm}, 1);
    testGeometry1x1x1<MeshFrameworkAudit, MeshFramework>(mesh);
  }
}


TEST(MESH_GEOMETRY_1BOX_EXO)
{
  // only works in MSTK or MOAB
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }
  if (framework_enabled(Framework::MOAB)) {
    frameworks.push_back(Framework::MOAB);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 1x1x1 Exo geometry with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createFrameworkUnstructured(Preference{frm}, "test/hex_1x1x1_sets.exo");
    testGeometry1x1x1<MeshFrameworkAudit, MeshFramework>(mesh);
  }
}


TEST(MESH_GEOMETRY_2BOX)
{
  // only works in MSTK or Simple
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }

  if (getDefaultComm()->NumProc() == 1)
    frameworks.push_back(Framework::SIMPLE);

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 2x2x2 with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createFrameworkStructuredUnitCube(Preference{frm}, 2);
    testGeometry2x2x2<MeshFrameworkAudit, MeshFramework>(mesh);
  }
}


TEST(MESH_GEOMETRY_3BOX_EXO)
{
  // only works in MSTK or MOAB
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }
  if (framework_enabled(Framework::MOAB)) {
    frameworks.push_back(Framework::MOAB);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Box 3x3x3 Exo with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createFrameworkUnstructured(Preference{frm}, "test/hex_3x3x3_sets.exo");
    testMeshAudit<MeshFrameworkAudit, MeshFramework>(mesh);
  }
}


TEST(MESH_GEOMETRY_FRACTURE_EXO)
{
  // only works in MSTK or MOAB
  std::vector<Framework> frameworks;
  if (framework_enabled(Framework::MSTK)) {
    frameworks.push_back(Framework::MSTK);
  }
  if (framework_enabled(Framework::MOAB)) {
    frameworks.push_back(Framework::MOAB);
  }

  for (const auto& frm : frameworks) {
    std::cout << std::endl
              << "Testing 3D Fracture Exo with " << AmanziMesh::framework_names.at(frm) << std::endl
              << "------------------------------------------------" << std::endl;
    auto mesh = createFrameworkUnstructured(Preference{frm}, "test/fractures.exo");
    testMeshAudit<MeshFrameworkAudit, MeshFramework>(mesh);
  }
}
