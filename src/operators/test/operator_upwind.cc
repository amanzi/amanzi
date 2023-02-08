/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "MeshFactory.hh"
#include "Mesh_Algorithms.hh"
#include "VerboseObject.hh"

// Operators
#include "OperatorDefs.hh"
#include "OperatorUtils.hh"
#include "UpwindDivK.hh"
#include "UpwindGravity.hh"
#include "UpwindFlux.hh"
#include "UpwindFluxAndGravity.hh"

namespace Amanzi {

class Model {
 public:
  Model(Teuchos::RCP<const AmanziMesh::Mesh> mesh) : mesh_(mesh){};
  ~Model(){};

  // main members
  double Value(int c, double pc) const { return analytic(pc); }

  double analytic(double pc) const { return 1e-5 + pc; }

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

} // namespace Amanzi


/* *****************************************************************
* Test one upwind model.
* **************************************************************** */
using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Operators;

template <class UpwindClass>
void
RunTestUpwind(std::string method)
{
  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: 1st-order convergence for upwind \"" << method << "\"\n";

  // read parameter list
  std::string xmlFileName = "test/operator_upwind.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  Teuchos::ParameterList plist = xmlreader.getParameters();

  // create an SIMPLE mesh framework
  Teuchos::ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  for (int n = 4; n < 17; n *= 2) {
    Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n, n);

    int ncells_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::ALL);
    int nfaces_owned = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
    int nfaces_wghost = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::ALL);

    // create model of nonlinearity
    Teuchos::RCP<Model> model = Teuchos::rcp(new Model(mesh));

    // create boundary data
    std::vector<int> bc_model(nfaces_wghost, OPERATOR_BC_NONE);
    std::vector<double> bc_value(nfaces_wghost);
    for (int f = 0; f < nfaces_wghost; f++) {
      const Point& xf = mesh->getFaceCentroid(f);
      if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6 ||
          fabs(xf[1] - 1.0) < 1e-6 || fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6)

        bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = model->analytic(xf[0]);
    }

    // create and initialize cell-based field
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);

    CompositeVector field(*cvs);
    Epetra_MultiVector& fcells = *field.ViewComponent("cell", true);
    Epetra_MultiVector& ffaces = *field.ViewComponent("face", true);

    for (int c = 0; c < ncells_wghost; c++) {
      const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
      fcells[0][c] = model->Value(c, xc[0]);
    }

    // add boundary face component
    for (int f = 0; f != bc_model.size(); ++f) {
      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);
        ffaces[0][f] = model->Value(c, bc_value[f]);
      }
    }

    // create and initialize face-based flux field
    cvs = CreateCompositeVectorSpace(mesh, "face", AmanziMesh::Entity_kind::FACE, 1, true);

    CompositeVector flux(*cvs), solution(*cvs);
    Epetra_MultiVector& u = *flux.ViewComponent("face", true);

    Point vel(1.0, 2.0, 3.0);
    for (int f = 0; f < nfaces_wghost; f++) {
      const Point& normal = mesh->getFaceNormal(f);
      u[0][f] = vel * normal;
    }

    // Create two upwind models
    Teuchos::ParameterList& ulist = plist.sublist("upwind");
    UpwindClass upwind(mesh);
    upwind.Init(ulist);

    upwind.Compute(flux, solution, bc_model, field);

    // calculate errors
    double error(0.0);
    for (int f = 0; f < nfaces_owned; f++) {
      const Point& xf = mesh->getFaceCentroid(f);
      double exact = model->analytic(xf[0]);

      error += pow(exact - ffaces[0][f], 2.0);

      CHECK(ffaces[0][f] >= 0.0);
    }
#ifdef HAVE_MPI
    double tmp = error;
    mesh->getComm()->SumAll(&tmp, &error, 1);
    int itmp = nfaces_owned;
    mesh->getComm()->SumAll(&itmp, &nfaces_owned, 1);
#endif
    error = sqrt(error / nfaces_owned);

    if (comm->MyPID() == 0) printf("n=%2d %s=%8.4f\n", n, method.c_str(), error);
  }
}

TEST(UPWIND_FLUX)
{
  RunTestUpwind<UpwindFlux>("flux");
}

TEST(UPWIND_DIVK)
{
  RunTestUpwind<UpwindDivK>("divk");
}

TEST(UPWIND_GRAVITY)
{
  RunTestUpwind<UpwindGravity>("gravity");
}

// TEST(UPWIND_FLUX_GRAVITY) {
//  RunTestUpwind<UpwindFluxAndGravity>("flux_gravity");
// }
