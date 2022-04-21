/*
  State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Tests for state as a container of data
  NOTE: this test passes if it compiles!
*/

// TPLs
#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

// Amanzi
#include "errors.hh"
#include "MeshFactory.hh"

// Amanzi::State
#include "IO.hh"
#include "Data.hh"
#include "Data_Helpers.hh"
#include "State.hh"

struct MyPoint {
  double a;
  double b;
};

using MyPointList = std::vector<MyPoint>;

bool inline UserInitialize(Teuchos::ParameterList& plist, MyPointList& t,
                           const Amanzi::Key& fieldname,
                           const std::vector<std::string>& subfieldnames) {
  std::cout << "found it!" << std::endl;
  return true;
}

void UserWriteVis(const Amanzi::Visualization& vis,
                  const Amanzi::Key& fieldname,
                  const std::vector<std::string>& subfieldnames,
                  const MyPointList& vec) {}

void UserWriteCheckpoint(const Amanzi::Checkpoint& chkp,
                         const Amanzi::Key& fieldname,
                         const std::vector<std::string>& subfieldnames,
                         const MyPointList& vec) {
}
void UserReadCheckpoint(const Amanzi::Checkpoint& chkp,
                        const Amanzi::Key& fieldname,
                        const std::vector<std::string>& subfieldnames,
                        MyPointList& vec) {}

TEST(STATE_EXTENSIBILITY_CREATION) {
  using namespace Amanzi;

  auto comm = Amanzi::getDefaultComm();
  Teuchos::ParameterList region_list;
  auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  Amanzi::AmanziMesh::Preference pref;
  pref.clear();
  pref.push_back(Amanzi::AmanziMesh::Framework::MSTK);

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> m = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1);

  std::string xmlFileName = "test/state_extensibility.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  auto plist = Teuchos::parameterList(xmlreader.getParameters());

  State s(*Teuchos::sublist(plist, "state"));
  s.RegisterDomainMesh(m);
  s.Require<MyPointList>("my_points", Tags::DEFAULT, "my_points");
  s.GetRecordW("my_points", "my_points").set_io_vis();
  s.Setup();
  s.InitializeFields();

  Visualization vis(plist->sublist("visualization"));
  vis.set_mesh(m);
  vis.CreateFiles();
  WriteVis(vis, s);

  Checkpoint chkp(plist->sublist("checkpoint"), s);
  Amanzi::WriteCheckpoint(chkp, comm, s, 0.0);
}

