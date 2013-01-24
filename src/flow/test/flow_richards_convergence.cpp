/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "Epetra_SerialComm.h"
#include "Epetra_MpiComm.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Richards_PK.hpp"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::AmanziFlow;


/* ******************************************************************
* Calculate L2 error in pressure.                                                    
****************************************************************** */
double calculatePressureCellError(Teuchos::RCP<Mesh> mesh, Epetra_Vector& pressure)
{
  double k1 = 0.5, k2 = 2.0, g = 2.0, a = 5.0, cr = 1.02160895462971866;  // analytical data
  double f1 = sqrt(1.0 - g * k1 / cr);
  double f2 = sqrt(g * k2 / cr - 1.0);

  double pressure_exact, error_L2 = 0.0;
  for (int c = 0; c < pressure.MyLength(); c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);

    double z = xc[2];
    if (z < -a)
      pressure_exact = f1 * tan(cr * (z + 2*a) * f1 / k1);
    else
      pressure_exact = -f2 * tanh(cr * f2 * (z + a) / k2 - atanh(f1 / f2 * tan(cr * a * f1 / k1)));
    // cout << z << " " << pressure[c] << " exact=" <<  pressure_exact << endl;
    error_L2 += std::pow(pressure[c] - pressure_exact, 2.0) * volume;
  }
  return sqrt(error_L2);
}


/* ******************************************************************
* Calculate l2 error (small l) in darcy flux.                                                    
****************************************************************** */
double calculateDarcyFluxError(Teuchos::RCP<Mesh> mesh, Epetra_Vector& darcy_flux)
{
  double cr = 1.02160895462971866;  // analytical data
  AmanziGeometry::Point velocity_exact(0.0, 0.0, -cr);

  int nfaces = darcy_flux.MyLength();
  double error_l2 = 0.0;
  for (int f = 0; f < nfaces; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    // cout << f << " " << darcy_flux[f] << " exact=" << velocity_exact * normal << endl;
    error_l2 += std::pow(darcy_flux[f] - velocity_exact * normal, 2.0);
  }
  return sqrt(error_l2 / nfaces);
}


/* ******************************************************************
* Calculate L2 divergence error in darcy flux.                                                    
****************************************************************** */
double calculateDarcyDivergenceError(Teuchos::RCP<Mesh> mesh, Epetra_Vector& darcy_flux)
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
      div += darcy_flux[f] * dirs[i];
    }
    error_L2 += div*div / mesh->cell_volume(c);
  }
  return sqrt(error_L2);
}


TEST(FLOW_RICHARDS_CONVERGENCE) {
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout <<"Richards: convergence Analysis" << std::endl;

  Teuchos::ParameterList parameter_list;
  string xmlFileName = "test/flow_richards_convergence.xml";
  // DEPRECATED  updateParametersFromXmlFile(xmlFileName, &parameter_list);
  
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  parameter_list = xmlreader.getParameters();

  // convergence estimate
  std::vector<double> h, p_error, v_error;

  for (int n = 40; n < 321; n*=2) {
    Teuchos::ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
    GeometricModelPtr gm = new GeometricModel(3, region_list, comm);
    
    FrameworkPreference pref;
    pref.clear();
    pref.push_back(MSTK);

    MeshFactory meshfactory(comm);
    meshfactory.preference(pref);
    Teuchos::RCP<Mesh> mesh = meshfactory(0.0, 0.0, -10.0, 1.0, 1.0, 0.0, 1, 1, n, gm);

    // create Richards process kernel
    Teuchos::ParameterList state_list = parameter_list.get<Teuchos::ParameterList>("State");
    State* S = new State(state_list, mesh);
    S->set_time(0.0);

    Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(*S));
    Richards_PK* RPK = new Richards_PK(parameter_list, FS);

    RPK->InitPK();  // setup the problem
    RPK->InitSteadyState(0.0, 0.01);
    if (n == 40) RPK->PrintStatistics();
    RPK->set_verbosity(FLOW_VERBOSITY_NONE);

    RPK->AdvanceToSteadyState(0.0, 0.01);
    RPK->CommitState(FS);

    double pressure_err, flux_err, div_err;  // error checks
    pressure_err = calculatePressureCellError(mesh, FS->ref_pressure());
    flux_err = calculateDarcyFluxError(mesh, FS->ref_darcy_flux());
    div_err = calculateDarcyDivergenceError(mesh, FS->ref_darcy_flux());

    if (n == 80) CHECK(pressure_err < 5.0e-2 && flux_err < 5.0e-2);
    int num_nonlinear_steps = -1;
    printf("n=%3d itrs=%4d  L2_pressure_err=%7.3e  l2_flux_err=%7.3e  L2_div_err=%7.3e\n",
        n, num_nonlinear_steps, pressure_err, flux_err, div_err);

    delete RPK;
    delete S;

    h.push_back(10.0 / n);
    p_error.push_back(pressure_err);
    v_error.push_back(flux_err);
  }

  // convergence rates
  double p_rate = Amanzi::AmanziFlow::bestLSfit(h, p_error);
  double v_rate = Amanzi::AmanziFlow::bestLSfit(h, v_error);
  printf("convergence rates: %23.2f %22.2f\n", p_rate, v_rate);

  CHECK_CLOSE(p_rate, 2.0, 0.2);
  CHECK_CLOSE(v_rate, 2.0, 0.2);

  delete comm;
}

