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

  MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(Preference({ Framework::MSTK, Framework::STK }));

  for (int n = 4; n < 17; n *= 2) {
    std::string setname("fractures");
    auto mesh3D = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, n, n, n, true, true);
    Teuchos::RCP<const Mesh> mesh = Teuchos::rcp(
      new MeshExtractedManifold(mesh3D, setname, AmanziMesh::FACE, comm, gm, plist, true, false));

    int ncells_wghost = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
    int nfaces_owned = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
    int nfaces_wghost = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

    // create boundary data
    std::vector<int> bc_model(nfaces_wghost, OPERATOR_BC_NONE);
    std::vector<double> bc_value(nfaces_wghost);
    for (int f = 0; f < nfaces_wghost; f++) {
      const Point& xf = mesh->face_centroid(f);
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
      const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
      field_c[0][c] = Value(xc);
    }

    // add boundary face component
    int ndir(0);
    const auto& fmap = *field.ComponentMap("face", true);
    for (int f = 0; f != nfaces_owned; ++f) {
      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        int g = fmap.FirstPointInElement(f);
        field_f[0][g] = bc_value[f];
        ndir++;
      }
    }

    // create and initialize face-based flux field
    AmanziMesh::Entity_ID_List cells;

    CompositeVector flux(*cvs), solution(*cvs);
    auto& flux_f = *flux.ViewComponent("face", true);

    int dir;
    Point vel(1.0, 2.0, 3.0);
    for (int f = 0; f < nfaces_owned; f++) {
      int g = fmap.FirstPointInElement(f);
      int ndofs = fmap.ElementSize(f);

      mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      AMANZI_ASSERT(cells.size() == ndofs);

      for (int i = 0; i < ndofs; ++i) {
        const Point& normal = mesh->face_normal(f, false, cells[i], &dir);
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
      mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

      const Point& xf = mesh->face_centroid(f);
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
    mesh->get_comm()->SumAll(&tmp, &error, 1);
    int itmp = ndir;
    mesh->get_comm()->SumAll(&itmp, &ndir, 1);
    itmp = nfaces_owned;
    mesh->get_comm()->SumAll(&itmp, &nfaces_owned, 1);
#endif
    error /= nfaces_owned;

    if (comm->MyPID() == 0) printf("n=%2d  dirichlet faces=%5d  error=%7.4f\n", n, ndir, error);
  }
}
