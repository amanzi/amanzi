/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/*
 State

 Tests for state as a container of data

 NOTE: this test passes if it compiles!
*/

// TPLs
#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"

#include "AmanziComm.hh"
#include "MeshFactory.hh"

#include "MeshFactory.hh"
#include "State.hh"
#include "errors.hh"

struct MyPoint {
  double a;
  double b;
};

using MyPointList = std::vector<MyPoint>;

bool inline UserInitialize(Teuchos::ParameterList& plist,
                           const Teuchos::ParameterList& attrs, MyPointList& t)
{
  std::cout << "Successfully initialized a MyPoint!" << std::endl;
  return true;
}

void
UserWriteVis(const Amanzi::Visualization& vis,
             const Teuchos::ParameterList& attrs, const MyPointList& vec)
{}
void
UserWriteCheckpoint(const Amanzi::Checkpoint& chkp,
                    const Teuchos::ParameterList& attrs, const MyPointList& vec)
{}
void
UserReadCheckpoint(const Amanzi::Checkpoint& chkp,
                   const Teuchos::ParameterList& attrs, MyPointList& vec)
{}

TEST(STATE_EXTENSIBILITY_CREATION)
{
  using namespace Amanzi;

  auto comm = Amanzi::getDefaultComm();
  Teuchos::ParameterList region_list;
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm = Teuchos::rcp(
    new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
  auto m = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1);

  std::string xmlFileName = "test/state_extensibility.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  auto plist = Teuchos::parameterList(xmlreader.getParameters());

  State s(Teuchos::sublist(plist, "state"));
  s.RegisterDomainMesh(m);
  s.Require<MyPointList>("my_points", "", "my_points");
  s.GetRecordW("my_points", "my_points").set_io_vis();
  s.Setup();

  Visualization vis(Teuchos::sublist(plist, "visualization"), m);
  WriteVis(vis, s);

  Checkpoint chkp(Teuchos::sublist(plist, "checkpoint"), comm);
  WriteCheckpoint(chkp, s);
}
