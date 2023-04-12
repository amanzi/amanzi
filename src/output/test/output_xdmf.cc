/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
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


SUITE(OUTPUT_XDMF) {
  TEST_FIXTURE(output_test_harness, XDMF_STRUCTURED)
  {
    create_mesh_structured();
    create_data();

    // Write a file which contains both mesh and data.
    Teuchos::ParameterList plist;
    std::string fnb = std::string("amanzi_vis_xdmf_structured_np") + std::to_string(comm->getSize());
    plist.set("file name base", fnb);
    Amanzi::OutputXDMF io(plist, mesh);

    test_write_multiple(io);

    if (comm->getRank() == 0) {
      std::string prefix = std::string("amanzi_vis_xdmf_structured_np") + std::to_string(comm->getSize());
      std::string test_prefix = std::string("test/gold_") + prefix;

      // mesh h5 file
      std::string command = std::string("h5diff ") + prefix+"_mesh.h5 test/gold_"+prefix+"_mesh.h5";
      int res = system(command.c_str());
      std::cout << "res = " << res << std::endl;
      CHECK(!res);

      // data h5 file
      command = std::string("h5diff ") + prefix+"_data.h5 test/gold_"+prefix+"_data.h5";
      res = system(command.c_str());
      std::cout << "res = " << res << std::endl;
      CHECK(!res);

      // xmf files
      CHECK(compareTextFiles(prefix+"_mesh.VisIt.xmf", test_prefix+"_mesh.VisIt.xmf"));
      CHECK(compareTextFiles(prefix+"_data.VisIt.xmf", test_prefix+"_data.VisIt.xmf"));
      CHECK(compareTextFiles(prefix+"_mesh.h5.0.xmf", test_prefix+"_mesh.h5.0.xmf"));
      CHECK(compareTextFiles(prefix+"_data.h5.0.xmf", test_prefix+"_data.h5.0.xmf"));
      CHECK(compareTextFiles(prefix+"_data.h5.3.xmf", test_prefix+"_data.h5.3.xmf"));
    }

  }


  TEST_FIXTURE(output_test_harness, XDMF_POLYHEDRAL)
  {
    create_mesh_polyhedral();
    create_data();

    // Write a file which contains both mesh and data.
    Teuchos::ParameterList plist;
    std::string fnb = std::string("amanzi_vis_xdmf_polyhedral_np") + std::to_string(comm->getSize());
    plist.set("file name base", fnb);
    Amanzi::OutputXDMF io(plist, mesh);

    test_write_multiple(io);

    if (comm->getRank() == 0) {
      std::string prefix = std::string("amanzi_vis_xdmf_polyhedral_np") + std::to_string(comm->getSize());
      std::string test_prefix = std::string("test/gold_") + prefix;

      // mesh h5 file
      std::string command = std::string("h5diff ") + prefix+"_mesh.h5 test/gold_"+prefix+"_mesh.h5";
      int res = system(command.c_str());
      std::cout << "res = " << res << std::endl;
      CHECK(!res);

      // data h5 file
      command = std::string("h5diff ") + prefix+"_data.h5 test/gold_"+prefix+"_data.h5";
      res = system(command.c_str());
      std::cout << "res = " << res << std::endl;
      CHECK(!res);

      // xmf files
      CHECK(compareTextFiles(prefix+"_mesh.VisIt.xmf", test_prefix+"_mesh.VisIt.xmf"));
      CHECK(compareTextFiles(prefix+"_data.VisIt.xmf", test_prefix+"_data.VisIt.xmf"));
      CHECK(compareTextFiles(prefix+"_mesh.h5.0.xmf", test_prefix+"_mesh.h5.0.xmf"));
      CHECK(compareTextFiles(prefix+"_data.h5.0.xmf", test_prefix+"_data.h5.0.xmf"));
      CHECK(compareTextFiles(prefix+"_data.h5.3.xmf", test_prefix+"_data.h5.3.xmf"));
    }
  }

}
