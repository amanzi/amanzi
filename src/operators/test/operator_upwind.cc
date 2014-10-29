/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "GMVMesh.hh"

#include "OperatorDefs.hh"
#include "Upwind.hh"

namespace Amanzi{

class Model {
 public:
  Model(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh) {};
  ~Model() {};

  // main members
  double Value(int c, double pc, std::string name) const { 
    return analytic(pc); 
  }

  double Value(int c, double pc) const { 
    return analytic(pc); 
  }

  double analytic(double pc) const { return 1e-5 + pc; }

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

typedef double(Model::*ModelUpwindFn)(int c, double pc) const; 

}  // namespace Amanzi


/* *****************************************************************
* This test replaces diffusion tensor and boundary conditions by
* continuous functions.
* **************************************************************** */
TEST(UPWIND) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int MyPID = comm.MyPID();

  if (MyPID == 0) std::cout << "\nTest: Upwind models" << std::endl;

  // read parameter list
  std::string xmlFileName = "test/operator_upwind.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, &comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);
  pref.push_back(STKMESH);

  MeshFactory meshfactory(&comm);
  meshfactory.preference(pref);

  for (int n = 4; n < 17; n *= 2) {
    RCP<const Mesh> mesh = meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n, n, gm);

    int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
    int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
    int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

    // create model of nonlinearity
    RCP<Model> model = rcp(new Model(mesh));

    // create boundary data
    std::vector<int> bc_model(nfaces_wghost, OPERATOR_BC_NONE);
    std::vector<double> bc_value(nfaces_wghost);
    for (int f = 0; f < nfaces_wghost; f++) {
      const Point& xf = mesh->face_centroid(f);
      if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
          fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
          fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6) 

      bc_model[f] = OPERATOR_BC_FACE_DIRICHLET;
      bc_value[f] = model->analytic(xf[0]);
    }

    // create and initialize cell-based field 
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh);
    cvs->SetGhosted(true);
    cvs->SetComponent("cell", AmanziMesh::CELL, 1);
    cvs->SetOwned(false);
    cvs->AddComponent("face", AmanziMesh::FACE, 1);

    CompositeVector field(*cvs);
    Epetra_MultiVector& fcells = *field.ViewComponent("cell", true);
    Epetra_MultiVector& ffaces = *field.ViewComponent("face", true);

    for (int c = 0; c < ncells_wghost; c++) {
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      fcells[0][c] = model->Value(c, xc[0], " "); 
    }
    for (int f = 0; f < nfaces_wghost; f++) {
      const AmanziGeometry::Point& xf = mesh->face_centroid(f);
      ffaces[0][f] = model->Value(0, xf[0], " "); 
    }

    // create and initialize face-based flux field
    cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh);
    cvs->SetGhosted(true);
    cvs->SetComponent("face", AmanziMesh::FACE, 1);

    CompositeVector flux(*cvs), upw_field(*cvs);
    Epetra_MultiVector& u = *flux.ViewComponent("face", true);
  
    Point vel(1.0, 2.0, 3.0);
    for (int f = 0; f < nfaces_wghost; f++) {
      const Point& normal = mesh->face_normal(f);
      u[0][f] = vel * normal;
    }

    // Create upwind
    Upwind<Model> upwind(mesh, model);

    ParameterList& ulist = plist.sublist("upwind");
    upwind.Init(ulist);
    // upwind.Compute(flux, bc_model, bc_value, field, upw_field, " ");
    ModelUpwindFn func = &Model::Value;
    upwind.Compute(flux, bc_model, bc_value, field, upw_field, func);

    // calculate error
    Epetra_MultiVector& upw = *upw_field.ViewComponent("face");
    double error(0.0);
    for (int f = 0; f < nfaces_owned; f++) {
      const Point& xf = mesh->face_centroid(f);
      double exact = model->analytic(xf[0]);
      error += pow(exact - upw[0][f], 2.0);
      CHECK(upw[0][f] >= 0.0);
    }
    error = sqrt(error / nfaces_owned);
  
    printf("error=%8.4f\n", error);
  }
}


