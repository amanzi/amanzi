#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"

#include "MeshAudit.hh"

#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"


TEST(MSTK_DEFORM_VOLS_2D)
{

  Teuchos::RCP<Epetra_MpiComm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));

  // Define a box region to capture bottom boundary
  Teuchos::ParameterList param_list;
  Teuchos::ParameterList& regions_list = param_list.sublist("Regions");
  Teuchos::ParameterList& botreg_list = regions_list.sublist("Bottom Region");
  Teuchos::ParameterList& botreg_def = botreg_list.sublist("Region: Box");
  Teuchos::Array<double> lo_coord = Teuchos::tuple(-1.6,-0.01);
  Teuchos::Array<double> hi_coord = Teuchos::tuple(1.6,0.01);
  botreg_def.set< Teuchos::Array<double> >("Low Coordinate",lo_coord);
  botreg_def.set< Teuchos::Array<double> >("High Coordinate",hi_coord);

  //  Teuchos::writeParameterListToXmlOStream(param_list,std::cout);

  if (comm->NumProc() > 1) return;

  // Create a geometric model from region spec

  Amanzi::AmanziGeometry::GeometricModelPtr gm = 
    new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, comm.get());

  // Generate a mesh consisting of 3x3 elements 

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> 
    mesh(new Amanzi::AmanziMesh::Mesh_MSTK(-1.5,0.0,1.5,1.0,5,5,comm.get(),gm));

  int nc = 
    mesh->num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);

  // Request target volume of 50% for bottom two cells in center column
  // The others are unconstrained except for a barrier of minimum volume

  std::vector<double> orig_volumes, target_volumes, target_weights, min_volumes;
  orig_volumes.reserve(nc);
  target_volumes.reserve(nc);
  target_weights.reserve(nc);
  min_volumes.reserve(nc);

  for (int i = 0; i < nc; i++) {    
    orig_volumes[i] = mesh->cell_volume(i);
    // target_volumes[i] = orig_volumes[i];
    target_volumes[i] = 0.0;
    min_volumes[i] = 0.25*orig_volumes[i];

    Amanzi::AmanziGeometry::Point ccen = mesh->cell_centroid(i);
    if (fabs(ccen[0]) < 1.0e-06)
      if (ccen[1] < 1.0)
        target_volumes[i] = 0.5*orig_volumes[i];
  }

  std::vector<std::string> fixed_setnames;
  fixed_setnames.push_back("Bottom Region");
  bool move_vertical = true;
  int status = mesh->deform(target_volumes,min_volumes,fixed_setnames,
                            move_vertical);
  CHECK(status);

  for (int i = 0; i < nc; i++) {
    if (target_volumes[i] > 0.0 && target_volumes[i] < orig_volumes[i]) {
      double voldiff = 
        (mesh->cell_volume(i)-target_volumes[i])/target_volumes[i];
    
      // Check if volume difference is with 5% of target volume
      CHECK_CLOSE(0,voldiff,5.0e-02);
    }

    // Check that we didn't fall below the minimum prescribed volume
    CHECK(mesh->cell_volume(i) > min_volumes[i]);
  }

  mesh->write_to_exodus_file("deformed2.exo");
}



TEST(MSTK_DEFORM_VOLS_3D)
{
  Teuchos::RCP<Epetra_MpiComm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));

  // Define a box region to capture bottom boundary
  Teuchos::ParameterList param_list;
  Teuchos::ParameterList& regions_list = param_list.sublist("Regions");
  Teuchos::ParameterList& botreg_list = regions_list.sublist("Bottom Region");
  Teuchos::ParameterList& botreg_def = botreg_list.sublist("Region: Box");
  Teuchos::Array<double> lo_coord = Teuchos::tuple(-1.6,-1.6,-0.01);
  Teuchos::Array<double> hi_coord = Teuchos::tuple(1.6,1.6,0.01);
  botreg_def.set< Teuchos::Array<double> >("Low Coordinate",lo_coord);
  botreg_def.set< Teuchos::Array<double> >("High Coordinate",hi_coord);

  //  Teuchos::writeParameterListToXmlOStream(param_list,std::cout);

  // Create a geometric model from region spec

  Amanzi::AmanziGeometry::GeometricModelPtr gm = 
    new Amanzi::AmanziGeometry::GeometricModel(3, regions_list, comm.get());

  // Generate a mesh consisting of 3x3 elements 

  std::string meshfilename = (comm->NumProc() == 1) ? "test/hex_5x5x5.exo" : "test/hex_5x5x5.par";
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> 
    mesh(new Amanzi::AmanziMesh::Mesh_MSTK(meshfilename.c_str(),comm.get(),gm));

  int nc = 
    mesh->num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::USED);

  // Request target volume of 50% for bottom two cells in center column
  // The others are unconstrained except for a barrier of minimum volume

  std::vector<double> orig_volumes, target_volumes, target_weights, min_volumes;
  orig_volumes.reserve(nc);
  target_volumes.reserve(nc);
  target_weights.reserve(nc);
  min_volumes.reserve(nc);

  for (int i = 0; i < nc; i++) {    
    orig_volumes[i] = mesh->cell_volume(i);
    target_volumes[i] = orig_volumes[i];
    // target_weights[i] = 0.25;
    target_volumes[i] = 0.0;
    min_volumes[i] = 0.5*orig_volumes[i];

    Amanzi::AmanziGeometry::Point ccen = mesh->cell_centroid(i);
    if (fabs(ccen[0]) < 1.0e-06 && fabs(ccen[1]) < 1.0e-06)
      if (ccen[2] < 1.0)
        target_volumes[i] = 0.5*orig_volumes[i];
  }

  std::vector<std::string> fixed_setnames;
  fixed_setnames.push_back("Bottom Region");
  bool move_vertical = true;
  int status = mesh->deform(target_volumes,min_volumes,fixed_setnames,
                            move_vertical);
  CHECK(status);

  for (int i = 0; i < nc; i++) {
    if (target_volumes[i] > 0.0 && target_volumes[i] < orig_volumes[i]) {
      double voldiff = (mesh->cell_volume(i)-target_volumes[i])/target_volumes[i];
    
      // Check if volume difference is with 5% of target volume
      CHECK_CLOSE(0,voldiff,5.0e-02);
    }

    // Check that we didn't fall below the minimum prescribed volume
    CHECK(mesh->cell_volume(i) > min_volumes[i]);
  }

  mesh->write_to_exodus_file("deformed3.exo");
}

