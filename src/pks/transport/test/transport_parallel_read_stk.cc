/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport

*/

#include "MeshFactory.hh"

#include "transport_parallel_read.hh"

TEST(ADVANCE_WITH_STK_PARALLEL_READ)
{
  runTest(Amanzi::AmanziMesh::Framework::STK);
}
