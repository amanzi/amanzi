/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include "stdlib.h"
#include "math.h"

// TPLs
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "hdf5.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "CycleDriver.hh"
#include "MeshAudit.hh"
#include "eos_reg.hh"
#include "Mesh.hh"
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "models_energy_reg.hh"
#include "models_flow_reg.hh"
#include "PK_Factory.hh"
#include "PK.hh"
#include "pks_energy_reg.hh"
#include "pks_flow_reg.hh"
#include "pks_mpc_reg.hh"
#include "pks_transport_reg.hh"
#include "State.hh"

void CreateApertureFile(int ncells, double time);

TEST(MPC_DRIVER_THERMAL_FLOW_MATRIX_FRACTURE_RICHARDS)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  Comm_ptr_type comm = Amanzi::getDefaultComm();

  // setup a piecewice linear solution with a jump
  std::string xmlInFileName = "test/mpc_coupled_thermal_flow_richards.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName);

  // For now create one geometric model from all the regions in the spec
  Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  // create mesh
  auto mesh_list = Teuchos::sublist(plist, "mesh", true);
  mesh_list->set<bool>("request edges", true);
  mesh_list->set<bool>("request faces", true);
  MeshFactory factory(comm, gm, mesh_list);
  factory.set_preference(Preference({ Framework::MSTK }));
  auto mesh = factory.create("test/single_fracture_tet.exo");

  // create dummy observation data object
  Amanzi::ObservationData obs_data;

  Teuchos::ParameterList state_plist = plist->sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));
  S->RegisterMesh("domain", mesh);

  // create additional mesh for fracture
  std::vector<std::string> names;
  names.push_back("fracture");
  auto mesh_fracture_framework = Teuchos::rcp(new MeshExtractedManifold(
    mesh, "fracture", AmanziMesh::Entity_kind::FACE, comm, gm, mesh_list));
  auto mesh_fracture = Teuchos::rcp(new Mesh(
    mesh_fracture_framework, Teuchos::rcp(new Amanzi::AmanziMesh::MeshAlgorithms()), mesh_list));

  S->RegisterMesh("fracture", mesh_fracture);

  // create aperture file on rank 0
  int ncells_tmp = mesh_fracture->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int ncells(ncells_tmp);
  comm->SumAll(&ncells_tmp, &ncells, 1);

  if (comm->MyPID() == 0) CreateApertureFile(ncells, 300.0);

  // run simulation
  Amanzi::CycleDriver cycle_driver(plist, S, comm, obs_data);
  cycle_driver.Go();
}


/* ******************************************************************
* Write data to an HDF5 file.
****************************************************************** */
void CreateApertureFile(int ncells, double time)
{
  hid_t hout = H5Fcreate("test/aperture_dynamic_test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t gout = H5Gcreate(hout, "fracture-aperture.cell.0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // times for aperture are hardcoded
  std::vector<double> aperture0(ncells, 4.6416e-4);
  std::vector<double> aperture1(ncells, 5.6416e-4);

  std::vector<double> times = { 0.0, time };

  hsize_t dims[2] = { (hsize_t)ncells, 1 };
  hid_t dataspace = H5Screate_simple(2, dims, NULL);
  hid_t dataset0 = H5Dcreate(gout, "0", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset0, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, aperture0.data());
  H5Dclose(dataset0);

  hid_t dataset1 = H5Dcreate(gout, "1", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, aperture1.data());
  H5Dclose(dataset1);

  hsize_t time_dims[1] = {2};
  hid_t time_dataspace = H5Screate_simple(1, time_dims, NULL);
  hid_t time_dataset = H5Dcreate(hout, "time", H5T_NATIVE_DOUBLE, time_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(time_dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, times.data());
  H5Dclose(time_dataset);

  H5Sclose(dataspace);
  H5Sclose(time_dataspace);
  H5Gclose(gout);
  H5Fclose(hout);
}

