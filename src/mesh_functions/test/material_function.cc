#include <map>
#include <iostream>
#include <string>
#include <vector>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "UnitTest++.h"

#include "errors.hh"
#include "ConstantFunction.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"

#include "MaterialMeshFunction.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Functions;

int main (int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return UnitTest::RunAllTests();
}


TEST(MESH2D)
{
  Epetra_MpiComm *comm;
  comm = new Epetra_MpiComm(MPI_COMM_WORLD);

  Teuchos::Array<double> corner_min(Teuchos::tuple(0.0, 0.0));
  Teuchos::Array<double> corner_max(Teuchos::tuple(0.2, 0.2));
  Teuchos::Array<double> normals(Teuchos::tuple(1.0, 0.0, 0.0, 1.0));

  Teuchos::ParameterList regions;
  regions.sublist("LEFT").sublist("region: box volume fractions")
      .set("corner coordinate", corner_min)
      .set("opposite corner coordinate", corner_max)
      .set("normals", normals);

  Teuchos::RCP<AmanziGeometry::GeometricModel> 
      gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(2, regions, comm));
  MeshFactory factory(comm);
  Teuchos::RCP<Mesh> mesh = factory(0.0, 0.0, 1.0, 1.0, 4, 4, gm);

  std::vector<std::string> rgns;
  rgns.push_back("LEFT");

  Teuchos::RCP<MultiFunction> f = Teuchos::rcp(new MultiFunction(Teuchos::rcp(new ConstantFunction(1.0))));
  Teuchos::RCP<MeshFunction::Domain> domain = Teuchos::rcp(new MeshFunction::Domain(rgns, AmanziMesh::CELL));
  Teuchos::RCP<MeshFunction::Spec> spec = Teuchos::rcp(new MeshFunction::Spec(domain, f));

  MaterialMeshFunction mmf(mesh);
  mmf.AddSpec(spec);

  delete comm;
}

