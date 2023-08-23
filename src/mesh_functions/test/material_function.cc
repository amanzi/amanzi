/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <map>
#include <iostream>
#include <string>
#include <vector>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "UnitTest++.h"

#include "errors.hh"
#include "FunctionConstant.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "VerboseObject_objs.hh"

#include "MaterialMeshFunction.hh"

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


// a test class to access data
class DomainFunction : public MaterialMeshFunction {
 public:
  DomainFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : MaterialMeshFunction(mesh){};
  ~DomainFunction(){};

  // access specs
  const std::map<AmanziMesh::Entity_ID, double> get_ids(const AmanziMesh::Entity_kind& kind)
  {
    std::vector<Teuchos::RCP<MaterialSpec>>& ms_list = *material_specs_[kind];
    Teuchos::RCP<MaterialSpec>& ms = *ms_list.begin();
    return *ms->second;
  }
};


TEST(MESH2D)
{
  auto comm = getDefaultComm();

  Teuchos::Array<double> corner_min(Teuchos::tuple(0.0, 0.0));
  Teuchos::Array<double> corner_max(Teuchos::tuple(0.3, 0.3));
  Teuchos::Array<double> normals(Teuchos::tuple(1.0, 0.0, 0.0, 1.0));

  Teuchos::ParameterList regions;
  regions.sublist("RGN1")
    .sublist("region: box volume fractions")
    .set("corner coordinate", corner_min)
    .set("opposite corner coordinate", corner_max)
    .set("normals", normals);
  regions.sublist("RGN2")
    .sublist("region: box")
    .set("low coordinate", corner_min)
    .set("high coordinate", corner_max);

  Teuchos::RCP<AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new AmanziGeometry::GeometricModel(2, regions, *comm));
  MeshFactory meshfactory(comm, gm);
  Teuchos::RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 4, 4);

  // test first region
  std::vector<std::string> rgns;
  rgns.push_back("RGN1");
  AmanziMesh::Entity_kind kind = AmanziMesh::CELL;

  Teuchos::RCP<MultiFunction> f =
    Teuchos::rcp(new MultiFunction(Teuchos::rcp(new FunctionConstant(1.0))));
  Teuchos::RCP<MeshFunction::Domain> domain = Teuchos::rcp(new MeshFunction::Domain(rgns, kind));
  Teuchos::RCP<MeshFunction::Spec> spec = Teuchos::rcp(new MeshFunction::Spec(domain, f));

  {
    DomainFunction df(mesh);
    df.AddSpec(spec);
    const std::map<AmanziMesh::Entity_ID, double>& ids = df.get_ids(kind);

    double vofs[4] = { 1.0, 0.2, 0.2, 0.04 };
    int n(0);
    for (std::map<AmanziMesh::Entity_ID, double>::const_iterator it = ids.begin(); it != ids.end();
         ++it) {
      CHECK_CLOSE(it->second, vofs[n++], 1e-10);
    }
  }

  // test scond region
  rgns.clear();
  rgns.push_back("RGN2");

  domain = Teuchos::rcp(new MeshFunction::Domain(rgns, kind));
  spec = Teuchos::rcp(new MeshFunction::Spec(domain, f));

  {
    DomainFunction df(mesh);
    df.AddSpec(spec);
    const std::map<AmanziMesh::Entity_ID, double>& ids = df.get_ids(kind);

    for (std::map<AmanziMesh::Entity_ID, double>::const_iterator it = ids.begin(); it != ids.end();
         ++it) {
      CHECK_CLOSE(it->second, 1.0, 1e-10);
    }
  }
}
