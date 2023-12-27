/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
 State

 This test is simply for developers to stuff their own data in and try to see
 if it can compile with State.
*/

// TPLs
#include "UnitTest++.h"

#include "errors.hh"
#include "MeshFactory.hh"
#include "Patch.hh"
#include "data/DataFactory.hh"

#include "State.hh"
#include "Data_Helpers.hh"


TEST(STATE_CONTAINS_MYDATA)
{
  using namespace Amanzi;

  Impl::DataFactory fac = Impl::dataFactory<MultiPatch<double>, MultiPatchSpace>();

  // // Create the geometric model
  // Teuchos::ParameterList regions;
  // std::vector<double> low = { 0., 0., 0. };
  // std::vector<double> high = { 4., 4., 4. };
  // regions.sublist("left")
  //   .sublist("region: box")
  //   .set<Teuchos::Array<double>>("low coordinate", low)
  //   .set<Teuchos::Array<double>>("high coordinate", {2., 4., 4.});
  // regions.sublist("right")
  //   .sublist("region: box")
  //   .set<Teuchos::Array<double>>("low coordinate", {2.,0.,0.} )
  //   .set<Teuchos::Array<double>>("high coordinate", high);
  // regions.sublist("point")
  //   .sublist("region: point")
  //   .set<Teuchos::Array<double>("coordinate", { 2., 2., 2. });
  // auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, regions, *comm));

  // MeshFactory meshfac(comm, gm);
  // auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

  // State s;
  // // s.RegisterDomainMesh(mesh);

  // // require data with factory
  // s.Require<Patch, PatchSpace>("my_patch", Tags::DEFAULT, "my_patch");
}
