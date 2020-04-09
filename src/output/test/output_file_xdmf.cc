/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

// Tests FileXDMF writer

#include "UnitTest++.h"

#include "AmanziComm.hh"
#include "FileXDMF.hh"
#include "output_test_utils.hh"

using namespace Amanzi;
TEST(XDMF_WRITER)
{
  FileXDMF xdmf("visdump", 1818, 800, 7200);
  xdmf.CreateTimestep(0.0, 0, true);
  xdmf.WriteField<double>("base_porosity", AmanziMesh::Entity_kind::CELL);
  xdmf.WriteFields<double>("darcy_velocity", 3, AmanziMesh::Entity_kind::CELL);
  xdmf.CloseTimestep(0.0, 0, true);

  xdmf.CreateTimestep(2.46406570841889092e-02, 104, false);
  xdmf.WriteField<double>("base_porosity", AmanziMesh::Entity_kind::CELL);
  xdmf.WriteFields<double>("darcy_velocity", 3, AmanziMesh::Entity_kind::CELL);
  xdmf.CloseTimestep(2.46406570841889092e-02, 104, false);

  std::cout << "CHECK: visdump_mesh.VisIt.xmf" << std::endl;
  std::string f1g = read_text_file("test/gold_visdump_mesh.VisIt.xmf");
  std::string f1 = read_text_file("visdump_mesh.VisIt.xmf");
  CHECK_EQUAL(f1g, f1);

  std::cout << "CHECK: visdump_mesh.h5.0.xmf" << std::endl;
  std::string f2g = read_text_file("test/gold_visdump_mesh.h5.0.xmf");
  std::string f2 = read_text_file("visdump_mesh.h5.0.xmf");
  CHECK_EQUAL(f2g, f2);

  std::cout << "CHECK: visdump_data.VisIt.xmf" << std::endl;
  std::string f3g = read_text_file("test/gold_visdump_data.VisIt.xmf");
  std::string f3 = read_text_file("visdump_data.VisIt.xmf");
  CHECK_EQUAL(f3g, f3);

  std::cout << "CHECK: visdump_data.h5.0.xmf" << std::endl;
  std::string f4g = read_text_file("test/gold_visdump_data.h5.0.xmf");
  std::string f4 = read_text_file("visdump_data.h5.0.xmf");
  CHECK_EQUAL(f4g, f4);

  std::cout << "CHECK: visdump_data.h5.104.xmf" << std::endl;
  std::string f5g = read_text_file("test/gold_visdump_data.h5.104.xmf");
  std::string f5 = read_text_file("visdump_data.h5.104.xmf");
  CHECK_EQUAL(f5g, f5);
}
