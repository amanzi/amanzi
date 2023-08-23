/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include <map>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MultiFunction.hh"
#include "FunctionConstant.hh"
#include "CompositeVectorFunction.hh"
#include "errors.hh"

#include "VerboseObject_objs.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Functions;

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();
  auto result = UnitTest::RunAllTests();
  Kokkos::finalize();
  return result;
}

struct another_reference_mesh {
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;
  Teuchos::RCP<GeometricModel> gm;
  std::string LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP, INVALID;

  another_reference_mesh()
  {
    LEFT = "LEFT";
    RIGHT = "RIGHT";
    FRONT = "FRONT";
    BACK = "BACK";
    BOTTOM = "BOTTOM";
    TOP = "TOP";
    INVALID = "INVALID";

    comm = getDefaultComm();
    // Brick domain corners and outward normals to sides
    Teuchos::Array<double> corner_min(Teuchos::tuple(0.0, 0.0, 0.0));
    Teuchos::Array<double> corner_max(Teuchos::tuple(4.0, 4.0, 4.0));
    Teuchos::Array<double> left(Teuchos::tuple(-1.0, 0.0, 0.0));
    Teuchos::Array<double> right(Teuchos::tuple(1.0, 0.0, 0.0));
    Teuchos::Array<double> front(Teuchos::tuple(0.0, -1.0, 0.0));
    Teuchos::Array<double> back(Teuchos::tuple(0.0, 1.0, 0.0));
    Teuchos::Array<double> bottom(Teuchos::tuple(0.0, 0.0, -1.0));
    Teuchos::Array<double> top(Teuchos::tuple(0.0, 0.0, 1.0));
    // Create the geometric model
    Teuchos::ParameterList regions;
    regions.sublist("LEFT").sublist("region: plane").set("point", corner_min).set("normal", left);
    regions.sublist("FRONT").sublist("region: plane").set("point", corner_min).set("normal", front);
    regions.sublist("BOTTOM")
      .sublist("region: plane")
      .set("point", corner_min)
      .set("normal", bottom);
    regions.sublist("RIGHT").sublist("region: plane").set("point", corner_max).set("normal", right);
    regions.sublist("BACK").sublist("region: plane").set("point", corner_max).set("normal", back);
    regions.sublist("TOP").sublist("region: plane").set("point", corner_max).set("normal", top);
    regions.sublist("DOMAIN")
      .sublist("region: box")
      .set("low coordinate", corner_min)
      .set("high coordinate", corner_max);

    gm = Teuchos::rcp(new GeometricModel(3, regions, *comm));
    // Create the mesh
    MeshFactory meshfactory(comm, gm);
    mesh = meshfactory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
  }
};


TEST_FIXTURE(another_reference_mesh, cv_function)
{
  // make the mesh function
  Teuchos::RCP<const Function> constant_func = Teuchos::rcp(new FunctionConstant(1.0));
  std::vector<Teuchos::RCP<const Function>> constant_funcs(1, constant_func);
  Teuchos::RCP<MultiFunction> vector_func = Teuchos::rcp(new MultiFunction(constant_funcs));

  std::vector<std::string> regions(1, "DOMAIN");
  Teuchos::RCP<MeshFunction::Domain> domainC =
    Teuchos::rcp(new MeshFunction::Domain(regions, AmanziMesh::Entity_kind::CELL));
  Teuchos::RCP<MeshFunction::Spec> specC =
    Teuchos::rcp(new MeshFunction::Spec(domainC, vector_func));

  // This type of support will eventually be added, at least to MSTK.
  Teuchos::RCP<MeshFunction::Domain> domainF =
    Teuchos::rcp(new MeshFunction::Domain(regions, AmanziMesh::Entity_kind::FACE));
  Teuchos::RCP<MeshFunction::Spec> specF =
    Teuchos::rcp(new MeshFunction::Spec(domainF, vector_func));

  Teuchos::RCP<MeshFunction> meshfunc = Teuchos::rcp(new MeshFunction(mesh));
  meshfunc->AddSpec(specC);
  meshfunc->AddSpec(specF);

  // couple the function to the location names
  std::vector<std::string> names(2);
  names[0] = "cell";
  names[1] = "face";
  CompositeVectorFunction cvfunc(meshfunc, names);

  // make the CV
  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::Entity_kind::CELL;
  locations[1] = AmanziMesh::Entity_kind::FACE;
  std::vector<int> num_dofs(2, 1);

  Teuchos::RCP<CompositeVectorSpace> cv_sp = Teuchos::rcp(new CompositeVectorSpace());
  cv_sp->SetMesh(mesh)->SetGhosted(false)->SetComponents(names, locations, num_dofs);
  Teuchos::RCP<CompositeVector> cv = Teuchos::rcp(new CompositeVector(*cv_sp));
  cv->PutScalar(0.0);

  // apply the function to the vector
  cvfunc.Compute(0.0, cv.ptr());

  // Check
  int ncells =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c != ncells; ++c) { CHECK_CLOSE(1.0, (*cv)("cell", 0, c), 0.0000001); }

  int nfaces =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f != nfaces; ++f) { CHECK_CLOSE(1.0, (*cv)("face", 0, f), 0.0000001); }
}
