/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

// Flow
#include "Richards_PK.hh"
#include "Richards_SteadyState.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Flow;

/* ******************************************************************
* Calculate L2 error in pressure.                                                    
****************************************************************** */
double calculatePressureCellError(Teuchos::RCP<const Mesh> mesh, const Epetra_MultiVector& p)
{
  double k1 = 0.5, k2 = 2.0, g = 2.0, a = 5.0, cr = 1.02160895462971866;  // analytical data
  double f1 = sqrt(1.0 - g * k1 / cr);
  double f2 = sqrt(g * k2 / cr - 1.0);

  double pexact, error_L2 = 0.0;
  for (int c = 0; c < p.MyLength(); c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);

    double z = xc[1];
    if (z < -a) {
      pexact = f1 * tan(cr * (z + 2*a) * f1 / k1);
    } else {
      pexact = -f2 * tanh(cr * f2 * (z + a) / k2 - atanh(f1 / f2 * tan(cr * a * f1 / k1)));
      // std::cout << z << " " << p[0][c] << " exact=" <<  pexact << std::endl;
    }
    error_L2 += std::pow(p[0][c] - pexact, 2.0) * volume;
  }
  return sqrt(error_L2);
}


/* ******************************************************************
* Calculate l2 error (small l) in darcy flux.                                                    
****************************************************************** */
double calculateDarcyFluxError(Teuchos::RCP<const Mesh> mesh, const Epetra_MultiVector& flux)
{
  double cr = 1.02160895462971866;  // analytical data
  AmanziGeometry::Point velocity_exact(0.0, -cr);

  int nfaces = flux.MyLength();
  double error_l2 = 0.0;
  for (int f = 0; f < nfaces; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    // std::cout << f << " " << flux[0][f] << " exact=" << velocity_exact * normal << std::endl;
    error_l2 += std::pow(flux[0][f] - velocity_exact * normal, 2.0);
  }
  return sqrt(error_l2 / nfaces);
}


/* ******************************************************************
* Calculate L2 divergence error in darcy flux.                                                    
****************************************************************** */
double calculateDarcyDivergenceError(Teuchos::RCP<const Mesh> mesh, const Epetra_MultiVector& flux)
{
  double error_L2 = 0.0;
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;

    mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    double div = 0.0;
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      div += flux[0][f] * dirs[i];
    }
    error_L2 += div*div / mesh->cell_volume(c);
    // std::cout << c << " div=" << div << " err=" << error_L2 << std::endl;
  }
  return sqrt(error_L2);
}


TEST(FLOW_RICHARDS_CONVERGENCE) {
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout <<"Convergence analysis on three random meshes" << std::endl;

  std::string xmlFileName = "test/flow_richards_random.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // convergence estimate
  int nmeshes = plist->get<int>("number of meshes", 1);
  std::vector<double> h, p_error, v_error;

  for (int n = 0; n < nmeshes; n++) {  // Use "n < 3" for the full test
    Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("Regions");
    GeometricModelPtr gm = new GeometricModel(2, region_list, comm);
    
    FrameworkPreference pref;
    pref.clear();
    pref.push_back(MSTK);
    pref.push_back(STKMESH);

    MeshFactory meshfactory(comm);
    meshfactory.preference(pref);
    // Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, -10.0, 1.0, 0.0, 5, 50, gm);
    Teuchos::RCP<const Mesh> mesh;
    if (n == 0) {
      mesh = meshfactory("test/random_mesh1.exo", gm);
    } else if (n == 1) {
      mesh = meshfactory("test/random_mesh2.exo", gm);
    } else if (n == 2) {
      mesh = meshfactory("test/random_mesh3.exo", gm);
    }

    /* create a simple state and populate it */
    Amanzi::VerboseObject::hide_line_prefix = false;

    Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("State");
    Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
    S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

    Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
    Richards_PK* RPK = new Richards_PK(plist, "Flow", S, soln);

    RPK->Setup();
    S->Setup();
    S->InitializeFields();
    S->InitializeEvaluators();

    // create Richards process kernel
    RPK->Initialize();
    S->CheckAllFieldsInitialized();

    // solver the problem
    TI_Specs ti_specs;
    ti_specs.T0 = 0.0;
    ti_specs.dT0 = 1.0;
    ti_specs.T1 = 1.0e+4;
    ti_specs.max_itrs = 400;

    AdvanceToSteadyState(S, *RPK, ti_specs, soln);
    RPK->CommitStep(0.0, 1.0);

    S->GetFieldData("darcy_flux")->ScatterMasterToGhosted("face");
    const Epetra_MultiVector& p = *S->GetFieldData("pressure")->ViewComponent("cell");
    const Epetra_MultiVector& flux = *S->GetFieldData("darcy_flux")->ViewComponent("face", true);

    double pressure_err, flux_err, div_err;  // error checks
    pressure_err = calculatePressureCellError(mesh, p);
    flux_err = calculateDarcyFluxError(mesh, flux);
    div_err = calculateDarcyDivergenceError(mesh, flux);

    int num_bdf1_steps = ti_specs.num_itrs;
    printf("mesh=%d bdf1_steps=%d  L2_pressure_err=%7.3e  l2_flux_err=%7.3e  L2_div_err=%7.3e\n",
        n, num_bdf1_steps, pressure_err, flux_err, div_err);

    CHECK(pressure_err < 2e-1 && flux_err < 2e-1 && div_err < 1e-7);

    GMV::open_data_file(*mesh, (std::string)"flow.gmv");
    GMV::start_data();
    GMV::write_cell_data(p, 0, "pressure");
    GMV::close_data_file();

    delete RPK;
  }

  delete comm;
}

