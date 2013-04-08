#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include <map>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "composite_function.hh"
#include "constant-function.hh"
#include "composite_vector_function.hh"
#include "errors.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Functions;

int main (int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  return UnitTest::RunAllTests ();
}

struct another_reference_mesh
{
  Epetra_MpiComm *comm;
  Teuchos::RCP<Mesh> mesh;
  GeometricModel *gm;
  std::string LEFT;
  std::string RIGHT;
  std::string FRONT;
  std::string BACK;
  std::string BOTTOM;
  std::string TOP;
  std::string INVALID;

  another_reference_mesh()
  {
    LEFT   = "LEFT";
    RIGHT  = "RIGHT";
    FRONT  = "FRONT";
    BACK   = "BACK";
    BOTTOM = "BOTTOM";
    TOP    = "TOP";
    INVALID = "INVALID";

    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
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
    regions.sublist("LEFT").sublist("Region: Plane").
        set("Location",corner_min).set("Direction",left);
    regions.sublist("FRONT").sublist("Region: Plane").
        set("Location",corner_min).set("Direction",front);
    regions.sublist("BOTTOM").sublist("Region: Plane").
        set("Location",corner_min).set("Direction",bottom);
    regions.sublist("RIGHT").sublist("Region: Plane").
        set("Location",corner_max).set("Direction",right);
    regions.sublist("BACK").sublist("Region: Plane").
        set("Location",corner_max).set("Direction",back);
    regions.sublist("TOP").sublist("Region: Plane").
        set("Location",corner_max).set("Direction",top);
    regions.sublist("DOMAIN").sublist("Region: Box").
      set("Low Coordinate", corner_min).set("High Coordinate", corner_max);

    gm = new GeometricModel(3,regions,comm);
    // Create the mesh
    MeshFactory mesh_fact(comm);
    mesh = mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2, gm);
  }
};

TEST_FIXTURE(another_reference_mesh, cv_function)
{

  // make the mesh function
  Teuchos::RCP<const Function> constant_func = Teuchos::rcp(new ConstantFunction(1.0));
  std::vector<Teuchos::RCP<const Function> > constant_funcs(1,constant_func);
  Teuchos::RCP<VectorFunction> vector_func =
    Teuchos::rcp(new CompositeFunction(constant_funcs));

  std::vector<std::string> regions(1, "DOMAIN");
  Teuchos::RCP<MeshFunction::Domain> domainC =
    Teuchos::rcp(new MeshFunction::Domain(regions, AmanziMesh::CELL));
  Teuchos::RCP<MeshFunction::Spec> specC =
    Teuchos::rcp(new MeshFunction::Spec(domainC, vector_func));

  // THIS OUGHT TO WORK but fails for strange caching reasons in STK mesh.
  // This type of support will eventually be added, at least to MSTK.
  Teuchos::RCP<MeshFunction::Domain> domainF =
    Teuchos::rcp(new MeshFunction::Domain(regions, AmanziMesh::FACE));
  Teuchos::RCP<MeshFunction::Spec> specF =
    Teuchos::rcp(new MeshFunction::Spec(domainF, vector_func));

  Teuchos::RCP<MeshFunction> meshfunc = Teuchos::rcp(new MeshFunction(mesh));
  meshfunc->AddSpec(specC);
  // meshfunc->AddSpec(specF);

  // couple the function to the location names
  // std::vector<std::string> names(2);
  // names[0] = "cell";
  // names[1] = "face";
  std::vector<std::string> names(1, "cell");
  CompositeVectorFunction cvfunc(meshfunc, names);

  // make the CV
  names.resize(2);
  names[0] = "cell";
  names[1] = "face";

  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::CELL;
  locations[1] = AmanziMesh::FACE;
  std::vector<int> num_dofs(2,1);

  Teuchos::RCP<CompositeVector> cv = Teuchos::rcp(new CompositeVector(mesh,
          names, locations, num_dofs, false));
  cv->CreateData();
  cv->PutScalar(0.0);

  // apply the function to the vector
  cvfunc.Compute(0.0, cv.ptr());

  // Check
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=ncells; ++c) {
    CHECK_CLOSE(1.0, (*cv)("cell", 0, c), 0.0000001);
  }

  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f!=nfaces; ++f) {
    CHECK_CLOSE(0.0, (*cv)("face", 0, f), 0.0000001);
  }
}
