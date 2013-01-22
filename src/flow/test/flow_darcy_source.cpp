/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "gmv_mesh.hh"

#include "State.hpp"
#include "Flow_State.hpp"
#include "Darcy_PK.hpp"


double K11(const Amanzi::AmanziGeometry::Point& x, double t) { 
  return (x[0] + 1) * (x[0] + 1) + x[1] * x[1];
}
double K22(const Amanzi::AmanziGeometry::Point& x, double t) { 
  return (x[0] + 1) * (x[0] + 1);
}
double K12(const Amanzi::AmanziGeometry::Point& x, double t) { 
  return -x[0] * x[1];
}
double pressure_exact(const Amanzi::AmanziGeometry::Point& x, double t) { 
  double xx = x[0];
  double yy = x[1];
  return xx * xx * xx * yy * yy + xx * sin(2 * xx * yy) * sin(2 * yy);
}
double source_exact(const Amanzi::AmanziGeometry::Point& x, double t) { 
  double xx = x[0];
  double yy = x[1];
  return xx * yy;
}
Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& x, double t) { 
  Amanzi::AmanziGeometry::Point v(x[0], x[1]);
  return v;
}

/* *****************************************************************
* This test replaves tensor and boundary conditions by continuous
* functions. It is most sensitive to error in location of problem
* coefficients.
* **************************************************************** */
TEST(FLOW_DARCY_SOURCE) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::AmanziFlow;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) cout << "Test: 2D steady-state full elliptic solver" << endl;

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/flow_darcy_source.xml";
  
  // DEPRECATED  updateParametersFromXmlFile(xmlFileName, &parameter_list);
  ParameterXMLFileReader xmlreader(xmlFileName);
  parameter_list = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(2, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);
  RCP<Mesh> mesh = meshfactory(0.0, 0.0, 1.0, 1.0, 18, 18, gm);
  // RCP<Mesh> mesh = meshfactory("test/median32x33.exo", gm);

  // create and populate fake flow state
  Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(mesh));
  FS->set_permeability(1.0, 1.0, "Material 1");
  FS->set_porosity(1.0);
  FS->set_fluid_viscosity(1.0);
  FS->set_fluid_density(1.0);
  FS->set_gravity(0.0);

  // create Darcy process kernel
  Darcy_PK* DPK = new Darcy_PK(parameter_list, FS);
  DPK->InitPK();
  DPK->InitSteadyState(0.0, 1e-8);

  // reset problem coefficients
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  std::vector<WhetStone::Tensor>& K = DPK->get_K();

  for (int c = 0; c < K.size(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    K[c](0, 0) = K11(xc, 0.0);
    K[c](1, 1) = K22(xc, 0.0);
    K[c](0, 1) = K12(xc, 0.0);
  }

  // update Dirichlet boundary data 
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  std::vector<bc_tuple>& bc_values = DPK->get_bc_values();

  for (int f = 0; f < nfaces; f++) {
    if (bc_values[f][0] > 0.0) {
      const Point& xf = mesh->face_centroid(f);
      bc_values[f][0] = pressure_exact(xf, 0.0);
    }
  }

  // create a problem
  DPK->AssembleMatrixMFD();

  // update right-hand side
  Teuchos::RCP<Epetra_Vector> rhs = DPK->get_matrix()->rhs();
  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);
    (*rhs)[c] += source_exact(xc, 0.0) * volume;
  }

  // steady-state solution
  DPK->AdvanceToSteadyState();
  DPK->CommitState(FS);

  // calculate errors
  Epetra_Vector& p = FS->ref_pressure();
  Epetra_Vector& flux = FS->ref_darcy_flux();

  double p_error = 0.0;
  for (int c = 0; c < ncells; c++) {
    const Point& xc = mesh->cell_centroid(c);
    double volume = mesh->cell_volume(c);
    p_error += std::pow(pressure_exact(xc, 0.0) - p[c], 2.0) * volume;
  }

  double flux_error = 0.0;
  for (int f = 0; f < nfaces; f++) {
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    const Point& xf = mesh->face_centroid(f);
    const AmanziGeometry::Point& velocity = velocity_exact(xf, 0.0);
    flux_error += std::pow(flux[f] - velocity * normal, 2.0);
  }
  cout << p_error << " " << flux_error << endl; 

  GMV::open_data_file(*mesh, (std::string)"flow.gmv");
  GMV::start_data();
  GMV::write_cell_data(p, 0, "pressure");
  GMV::close_data_file();

  delete DPK;
}
