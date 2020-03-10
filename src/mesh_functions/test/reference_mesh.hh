#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AmanziComm.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"


using namespace Amanzi;

struct reference_mesh {
  Comm_ptr_type comm;
  Teuchos::RCP<AmanziMesh::Mesh> mesh;

  reference_mesh()
  {
    comm = getDefaultComm();

    // Brick domain corners and outward normals to sides
    Teuchos::Array<double> corner_min(Teuchos::tuple(0.0, 0.0, 0.0));
    Teuchos::Array<double> corner_max(Teuchos::tuple(4.0, 4.0, 4.0));
    Teuchos::Array<double> top_corner_half(Teuchos::tuple(2.0, 4.0, 4.0));
    Teuchos::Array<double> bottom_corner_half(Teuchos::tuple(2.0, 0.0, 0.0));

    Teuchos::Array<double> left(Teuchos::tuple(-1.0, 0.0, 0.0));
    Teuchos::Array<double> right(Teuchos::tuple(1.0, 0.0, 0.0));
    Teuchos::Array<double> front(Teuchos::tuple(0.0, -1.0, 0.0));
    Teuchos::Array<double> back(Teuchos::tuple(0.0, 1.0, 0.0));
    Teuchos::Array<double> bottom(Teuchos::tuple(0.0, 0.0, -1.0));
    Teuchos::Array<double> top(Teuchos::tuple(0.0, 0.0, 1.0));

    // Create the geometric model
    Teuchos::ParameterList regions;
    regions.sublist("DOMAIN")
        .sublist("region: all");
    regions.sublist("LEFT")
      .sublist("region: plane")
      .set("point", corner_min)
      .set("normal", left);
    regions.sublist("FRONT")
      .sublist("region: plane")
      .set("point", corner_min)
      .set("normal", front);
    regions.sublist("BOTTOM")
      .sublist("region: plane")
      .set("point", corner_min)
      .set("normal", bottom);
    regions.sublist("RIGHT")
      .sublist("region: plane")
      .set("point", corner_max)
      .set("normal", right);
    regions.sublist("BACK")
      .sublist("region: plane")
      .set("point", corner_max)
      .set("normal", back);
    regions.sublist("TOP")
      .sublist("region: plane")
      .set("point", corner_max)
      .set("normal", top);
    regions.sublist("HALF1")
        .sublist("region: box")
        .set("low coordinate", corner_min)
        .set("high coordinate", top_corner_half);
    regions.sublist("HALF2")
        .sublist("region: box")
        .set("low coordinate", bottom_corner_half)
        .set("high coordinate", corner_max);
    Teuchos::RCP<AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new AmanziGeometry::GeometricModel(3, regions, *comm));
    // Create the mesh
    AmanziMesh::MeshFactory mesh_fact(comm, gm);
    mesh = mesh_fact.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
  }
};
