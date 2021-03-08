#include <UnitTest++.h>
#include <fstream>

#include "../Mesh_MSTK.hh"

#include "Epetra_Map.h"
#include "AmanziComm.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"

/*
 * Verification test for region inheritance from a volume-based body
 * onto a surface-based body.
 */
TEST(RegionInheritanceMSTK)
{

  std::string filename("test/hex_3x3x3_sets.exo");
  auto comm = Amanzi::getDefaultComm();

  std::vector<std::string> setnames;
  setnames.push_back(std::string("Top Surface"));
  setnames.push_back(std::string("Side Surface"));
  setnames.push_back(std::string("Union"));

  // Create param list in code
  Teuchos::ParameterList parameterlist;
  Teuchos::ParameterList& reg_spec = parameterlist.sublist("regions"); 
  
  // Define faceset 106 (top surface)
  Teuchos::ParameterList& top_surface = reg_spec.sublist(setnames[0]);
  Teuchos::ParameterList& top_surface_def = top_surface.sublist("region: labeled set");
  top_surface_def.set<std::string>("label","106");
  top_surface_def.set<std::string>("file",filename.c_str());
  top_surface_def.set<std::string>("format","Exodus II");
  top_surface_def.set<std::string>("entity","face");

  // Define faceset 102 (western face)
  Teuchos::ParameterList& side_surface = reg_spec.sublist(setnames[1]);
  Teuchos::ParameterList& side_surface_def = side_surface.sublist("region: labeled set");
  side_surface_def.set<std::string>("label","102");
  side_surface_def.set<std::string>("file",filename.c_str());
  side_surface_def.set<std::string>("format","Exodus II");
  side_surface_def.set<std::string>("entity","face");

  // Define union of regions
  Teuchos::Array<std::string> arrayregions;
  //std::string _label = Amanzi::AmanziInput::underscore("Top Surface");
  arrayregions.push_back(setnames[0]);
  arrayregions.push_back(setnames[1]);

  Teuchos::ParameterList& union_region = reg_spec.sublist(setnames[2]);
  Teuchos::ParameterList& union_region_def = union_region.sublist("region: logical");
  union_region_def.set<std::string>("operation","union");
  union_region_def.set("regions", arrayregions);

  // Initialize geometric model
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm;
  gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_spec, *comm));

  // Read the Exodus mesh from file
  auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::Mesh_MSTK(filename.c_str(),comm,gm));

  // Capture the cell IDs for the defined surface regions
  Amanzi::AmanziMesh::Entity_ID_List top_faces;
  Amanzi::AmanziMesh::Entity_ID_List side_faces;

  mesh->get_set_entities(setnames[0], Amanzi::AmanziMesh::Entity_kind::FACE, Amanzi::AmanziMesh::Parallel_type::ALL, &top_faces);
  mesh->get_set_entities(setnames[1], Amanzi::AmanziMesh::Entity_kind::FACE, Amanzi::AmanziMesh::Parallel_type::ALL, &side_faces);

  // Verify union of regions
  Amanzi::AmanziMesh::Entity_ID_List union_faces;
  mesh->get_set_entities(setnames[2], Amanzi::AmanziMesh::Entity_kind::FACE, Amanzi::AmanziMesh::Parallel_type::ALL, &union_faces);
  CHECK_EQUAL(top_faces.size() + side_faces.size(), union_faces.size());

  // Extract surface mesh
  Amanzi::AmanziMesh::Mesh_MSTK surfmesh(mesh,union_faces,Amanzi::AmanziMesh::Entity_kind::FACE, false, mesh->get_comm());

  Amanzi::AmanziMesh::Entity_ID_List surf_union_faces;
  surfmesh.get_set_entities(setnames[2], Amanzi::AmanziMesh::Entity_kind::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED, &surf_union_faces);

  // Verify that we have extracted the correct # of set entities
  CHECK_EQUAL(surf_union_faces.size(), union_faces.size());

  // Verify that surface mesh cell centroids == vol. mesh face centroids
  for (int i = 0; i < surf_union_faces.size(); i++) {
    Amanzi::AmanziGeometry::Point centroid1, centroid2;

    centroid1 = mesh->getFaceCentroid(union_faces[i])
    centroid2 = surfmesh.getCellCentroid(surf_union_faces[i])

    CHECK_EQUAL(abs(centroid1.x() - centroid2.x()) < 0.01, true);
    CHECK_EQUAL(abs(centroid1.y() - centroid2.y()) < 0.01, true);
    CHECK_EQUAL(abs(centroid1.z() - centroid2.z()) < 0.01, true);
  }

}



