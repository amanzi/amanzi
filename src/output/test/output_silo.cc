/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "UnitTest++.h"

#define MSTK_HAVE_MPI
#include "Mesh_MSTK.hh"

#include "../OutputSilo.hh"
#include "output_test_utils.hh"

SUITE(OUTPUT_SILO) {

  TEST_FIXTURE(output_test_harness, SILO_STRUCTURED)
  {
    create_mesh_structured();
    create_data();

    // Write a file which contains both mesh and data.
    Teuchos::ParameterList plist;
    std::string fnb = std::string("amanzi_vis_silo_structured_np") + std::to_string(comm->getSize());
    plist.set("file name base", fnb);
    Amanzi::OutputSilo io(plist, mesh);
    test_write_multiple(io);
  }


  TEST_FIXTURE(output_test_harness, SILO_POLYHEDRAL)
  {
    create_mesh_polyhedral();
    create_data();

    // Write a file which contains both mesh and data.
    Teuchos::ParameterList plist;
    std::string fnb = std::string("amanzi_vis_silo_polyhedral_np") + std::to_string(comm->getSize());
    plist.set("file name base", fnb);
    Amanzi::OutputSilo io(plist, mesh);

    test_write_multiple(io);
  }

}
