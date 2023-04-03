/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "MeshExtractedManifold.hh"
#include "MeshFactory.hh"
#include "Mesh_Algorithms.hh"
#include "VerboseObject.hh"

// Operators
#include "OperatorDefs.hh"
#include "OperatorUtils.hh"
#include "PDE_HelperDiscretization.hh"
#include "UpwindFluxManifolds.hh"

double
Value(const Amanzi::AmanziGeometry::Point& xyz)
{
  return 1e-5 + xyz[0] + 0 * xyz[1] + 0 * xyz[2];
}


/* *****************************************************************
* Test upwind method on manifolds
* **************************************************************** */
using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Operators;

TEST(UPWIND_FLUX_MANIFOLDS)
{
  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0)
    std::cout << "\nTest: 1st-order convergence for upwind \""
              << "\"\n";

  // read parameter list
  std::string xmlFileName = "test/operator_upwind.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // create an SIMPLE mesh framework
  Teuchos::ParameterList region_list = plist->sublist("regions");
  Teuchos::RCP<GeometricModel> gm = Teuchos::rcp(new GeometricModel(3, region_list, *comm));

  auto flist = Teuchos::rcp(new Teuchos::ParameterList());
  flist->set<bool>("request faces", true);  
  flist->set<bool>("request edges", true); 
  MeshFactory meshfactory(comm, gm, flist);
  meshfactory.set_preference(Preference({ Framework::MSTK }));

  for (int n = 4; n < 17; n *= 2) {
    std::string setname("fractures");
    auto mesh3D = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n, n);
    auto mesh_fw = Teuchos::rcp(
      new MeshExtractedManifold(mesh3D, setname, AmanziMesh::FACE, comm, gm, plist));
    auto mesh = Teuchos::rcp(new Mesh(mesh_fw, Teuchos::rcp(new MeshFrameworkAlgorithms()), Teuchos::null)); 

    int ncells_owned = mesh->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::OWNED);
    int ncells_wghost = mesh->getNumEntities(AmanziMesh::CELL, AmanziMesh::Parallel_kind::ALL);
    int nfaces_owned = mesh->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::OWNED);
    int nfaces_wghost = mesh->getNumEntities(AmanziMesh::FACE, AmanziMesh::Parallel_kind::ALL);

    // create boundary data
    std::vector<int> bc_model(nfaces_wghost, OPERATOR_BC_NONE);
    std::vector<double> bc_value(nfaces_wghost);
    for (int f = 0; f < nfaces_wghost; f++) {
      const Point& xf = mesh->getFaceCentroid(f);
      if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6 ||
          fabs(xf[1] - 1.0) < 1e-6 || fabs(xf[2]) < 1e-6 || fabs(xf[2] - 1.0) < 1e-6)
        bc_model[f] = OPERATOR_BC_DIRICHLET;
      bc_value[f] = Value(xf);
    }

    // create and initialize cell-based field
    auto cvs = Operators::CreateManifoldCVS(mesh);
    cvs->AddComponent("cell", AmanziMesh::CELL, 1);

    CompositeVector field(*cvs);
    auto& field_c = *field.ViewComponent("cell", true);
    auto& field_f = *field.ViewComponent("face", true);

    for (int c = 0; c < ncells_wghost; c++) {
      const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
      field_c[0][c] = Value(xc);
    }

    // add boundary face component
    int ndir(0);
    const auto& fmap = *field.ComponentMap("face", true);
    for (int f = 0; f != nfaces_owned; ++f) {
      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        int g = fmap.FirstPointInElement(f);
        int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);
        field_f[0][g] = bc_value[f];
        ndir++;
      }
    }

    // create and initialize face-based flux field

    CompositeVector flux(*cvs), solution(*cvs);
    auto& flux_f = *flux.ViewComponent("face", true);

    int dir;
    Point vel(1.0, 2.0, 3.0);
    for (int f = 0; f < nfaces_owned; f++) {
      int g = fmap.FirstPointInElement(f);
      int ndofs = fmap.ElementSize(f);

      auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == ndofs);

      for (int i = 0; i < ndofs; ++i) {
        const Point& normal = mesh->getFaceNormal(f, cells[i], &dir);
        flux_f[0][g + i] = vel * normal;
      }
    }

    // Create two upwind models
    Teuchos::ParameterList& ulist = plist->sublist("upwind");
    UpwindFluxManifolds upwind(mesh);
    upwind.Init(ulist);

    upwind.Compute(flux, bc_model, field);

    // calculate errors
    double error(0.0);
    for (int f = 0; f < nfaces_owned; f++) {
      auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);

      const Point& xf = mesh->getFaceCentroid(f);
      double exact = Value(xf);

      int g = fmap.FirstPointInElement(f);
      int ndofs = fmap.ElementSize(f);
      for (int i = 0; i < ndofs; ++i) {
        error += fabs(exact - field_f[0][g + i]);
        CHECK(field_f[0][g + i] > 0.0);
      }
    }
#ifdef HAVE_MPI
    double tmp = error;
    mesh->getComm()->SumAll(&tmp, &error, 1);
    int itmp = ndir;
    mesh->getComm()->SumAll(&itmp, &ndir, 1);
    itmp = nfaces_owned;
    mesh->getComm()->SumAll(&itmp, &nfaces_owned, 1);
#endif
    error /= nfaces_owned;

    if (comm->MyPID() == 0) printf("n=%2d  dirichlet faces=%5d  error=%7.4f\n", n, ndir, error);
  }
}
