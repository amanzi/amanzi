/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

// Tests OutputXDMF writer

#include "UnitTest++.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "output_test_utils.hh"

#include "AmanziComm.hh"
#include "AmanziVector.hh"
#include "OutputXDMF.hh"
#include "Mesh_MSTK.hh"

using namespace Amanzi;

TEST(WRITE_MESH)
{
  // current checks binary equality
  auto comm = Amanzi::getDefaultComm();
  auto mesh =
    Teuchos::rcp(new AmanziMesh::Mesh_MSTK("test/hillslope_noduff.exo", comm));

  Teuchos::ParameterList plist;
  plist.set("file name base", "visdump");

  auto poro = Teuchos::rcp(new Vector_type(mesh->cell_map(false)));
  poro->putScalar(0.5);

  auto darcy_velocity =
    Teuchos::rcp(new MultiVector_type(mesh->cell_map(false), 3));
  darcy_velocity->putScalar(1.0);

  {
    OutputXDMF out(plist, mesh);
    out.CreateFile(0.0, 0);

    {
      Teuchos::ParameterList plist("base_porosity");
      plist.set("location", AmanziMesh::CELL);
      out.Write(plist, *poro);
    }
    {
      Teuchos::ParameterList plist("darcy_velocity");
      plist.set("location", AmanziMesh::CELL);
      out.Write(plist, *darcy_velocity);
    }
    out.FinalizeFile();

    out.CreateFile(2.46406570841889092e-02, 104);
    {
      Teuchos::ParameterList plist("base_porosity");
      plist.set("location", AmanziMesh::CELL);
      out.Write(plist, *poro);
    }
    {
      Teuchos::ParameterList plist("darcy_velocity");
      plist.set("location", AmanziMesh::CELL);
      out.Write(plist, *darcy_velocity);
    }
    out.FinalizeFile();
  }

  if (comm->getSize() == 1) {
    int res = system("h5diff visdump_mesh.h5 test/gold_visdump_mesh.h5");
    std::cout << "res = " << res << std::endl;
    CHECK(!res);
  }
  if (comm->getRank() == 0) {
    CHECK(compareTextFiles("visdump_mesh.VisIt.xmf",
                           "test/gold_visdump_mesh.VisIt.xmf"));
    CHECK(compareTextFiles("visdump_data.VisIt.xmf",
                           "test/gold_visdump_data.VisIt.xmf"));
    CHECK(compareTextFiles("visdump_mesh.h5.0.xmf",
                           "test/gold_visdump_mesh.h5.0.xmf"));
    CHECK(compareTextFiles("visdump_data.h5.0.xmf",
                           "test/gold_visdump_data.h5.0.xmf"));
    CHECK(compareTextFiles("visdump_data.h5.104.xmf",
                           "test/gold_visdump_data.h5.104.xmf"));
  }
}
