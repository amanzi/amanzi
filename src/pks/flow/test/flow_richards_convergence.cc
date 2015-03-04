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

#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Richards_PK.hh"

#include "Richards_SteadyState.hh"


using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Flow;


/* ******************************************************************
* Calculate L2 error in pressure.                                                    
****************************************************************** */
double calculatePressureCellError(Teuchos::RCP<const Mesh> mesh,
                                  const Epetra_MultiVector& p)
{
  double k1 = 0.5, k2 = 2.0, g = 2.0, a = 5.0, cr = 1.02160895462971866;  // analytical data
  double f1 = sqrt(1.0 - g * k1 / cr);
  double f2 = sqrt(g * k2 / cr - 1.0);

  double pexact, error_L2 = 0.0;
  for (int c = 0; c < p.MyLength(); c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);

    double z = xc[2];
    if (z < -a) {
      pexact = f1 * tan(cr * (z + 2*a) * f1 / k1);
    } else {
      pexact = -f2 * tanh(cr * f2 * (z + a) / k2 - atanh(f1 / f2 * tan(cr * a * f1 / k1)));
    }
    error_L2 += std::pow(p[0][c] - pexact, 2.0) * volume;
    // std::cout << z << " " << p[0][c] << " exact=" <<  pexact << std::endl;
  }
  return sqrt(error_L2);
}


/* ******************************************************************
* Calculate l2 error (small l) in darcy flux.                                                    
****************************************************************** */
double calculateDarcyFluxError(Teuchos::RCP<const Mesh> mesh, 
                               const Epetra_MultiVector& flux)
{
  double cr = 1.02160895462971866;  // analytical data
  AmanziGeometry::Point velocity_exact(0.0, 0.0, -cr);

  int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  double error_l2 = 0.0;
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    error_l2 += std::pow(flux[0][f] - velocity_exact * normal, 2.0);
    // std::cout << f << " " << flux[0][f] << " exact=" << velocity_exact * normal << std::endl;
  }
  return sqrt(error_l2 / nfaces_owned);
}


/* ******************************************************************
* Calculate L2 divergence error in darcy flux.                                                    
****************************************************************** */
double calculateDarcyDivergenceError(Teuchos::RCP<const Mesh> mesh,
                                     const Epetra_MultiVector& flux)
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
    // std::cout << c << " div=" << div << " exact=0.0" << std::endl;
  }
  return sqrt(error_L2);
}


TEST(FLOW_RICHARDS_CONVERGENCE) {
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_WORLD);
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout <<"Richards: convergence Analysis" << std::endl;

  std::string xmlFileName = "test/flow_richards_convergence.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // convergence estimate
  std::vector<double> h, p_error, v_error;

  for (int n = 40; n < 161; n*=2) {
    Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("Regions");
    GeometricModelPtr gm = new GeometricModel(3, region_list, comm);
    
    FrameworkPreference pref;
    pref.clear();
    pref.push_back(MSTK);
    pref.push_back(STKMESH);

    MeshFactory meshfactory(comm);
    meshfactory.preference(pref);
    Teuchos::RCP<const Mesh> mesh = meshfactory(0.0, 0.0, -10.0, 1.0, 1.0, 0.0, 1, 1, n, gm);

    /* create and populate flow state */
    Amanzi::VerboseObject::hide_line_prefix = false;

    Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("State");
    Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
    S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

    /* create Richards process kernel */
    Richards_PK* RPK = new Richards_PK(plist, "Flow", S);
    RPK->Setup();
    S->Setup();
    S->InitializeFields();
    S->InitializeEvaluators();

    RPK->Initialize();
    S->CheckAllFieldsInitialized();
    RPK->InitTimeInterval();

    // solve the problem
    TI_Specs ti_specs;
    ti_specs.T0 = 0.0;
    ti_specs.dT0 = 1.0;
    ti_specs.T1 = 1.0e+4;
    ti_specs.max_itrs = 1000;

    AdvanceToSteadyState(S, *RPK, ti_specs, S->GetFieldData("pressure", "flow"));
    RPK->CommitState(0.0, S.ptr());

    double pressure_err, flux_err, div_err;  // error checks
    const Epetra_MultiVector& p = *S->GetFieldData("pressure")->ViewComponent("cell");
    const Epetra_MultiVector& flux = *S->GetFieldData("darcy_flux")->ViewComponent("face", true);

    pressure_err = calculatePressureCellError(mesh, p);
    flux_err = calculateDarcyFluxError(mesh, flux);
    div_err = calculateDarcyDivergenceError(mesh, flux);

    if (n == 80) CHECK(pressure_err < 5.0e-2 && flux_err < 5.0e-2);
    int num_nonlinear_steps = -1;
    printf("n=%3d itrs=%4d  L2_pressure_err=%7.3e  l2_flux_err=%7.3e  L2_div_err=%7.3e\n",
        n, num_nonlinear_steps, pressure_err, flux_err, div_err);

    delete RPK;

    h.push_back(10.0 / n);
    p_error.push_back(pressure_err);
    v_error.push_back(flux_err);
  }

  /* convergence rates */
  double p_rate = Amanzi::Flow::bestLSfit(h, p_error);
  double v_rate = Amanzi::Flow::bestLSfit(h, v_error);
  printf("convergence rates: %23.2f %22.2f\n", p_rate, v_rate);

  CHECK_CLOSE(p_rate, 2.0, 0.2);
  CHECK(v_rate > 1.8);

  delete comm;
}

