/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//! This is an example of using CompositeVectorFunctions to generate BCs

#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include <map>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MultiFunction.hh"
#include "CompositeVectorFunction.hh"
#include "errors.hh"
#include "CompositeVectorSpace.hh"
#include "DataStructuresHelpers.hh"

#include "test/reference_mesh.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Functions;

//
// This fixture is developer documentation, showing how to create boundary
// conditions from MultiPatches and MeshFunctions.
//
TEST_FIXTURE(reference_mesh, MESH_FUNCTION)
{
  std::string xmlFileName = "test/mesh_functions_bcs.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  auto plist = xmlreader.getParameters();

  // first, create a MultiPatch for all Dirichlet BCs
  int OPERATOR_BC_DIRICHLET = 1;
  MeshFunction mf_d(plist.sublist("dirichlet"),
                    mesh,
                    "boundary pressure",
                    Entity_kind::FACE,
                    OPERATOR_BC_DIRICHLET);
  auto mps = mf_d.createMPS(false);
  MultiPatch<double> dirichlet(mf_d.createMPS(false));
  mf_d.Compute(1.0, dirichlet);
  CHECK_EQUAL(OPERATOR_BC_DIRICHLET, (*dirichlet.space)[0]->flag_type);
  CHECK_EQUAL(OPERATOR_BC_DIRICHLET, dirichlet[0].space->flag_type);

  // next, one for Neumann
  int OPERATOR_BC_NEUMANN = 2;
  MeshFunction mf_n(
    plist.sublist("neumann"), mesh, "outward normal flux", Entity_kind::FACE, OPERATOR_BC_NEUMANN);
  MultiPatch<double> neumann(mf_n.createMPS(false));
  mf_n.Compute(1.0, neumann);
  CHECK_EQUAL(OPERATOR_BC_NEUMANN, (*neumann.space)[0]->flag_type);

  // next, one for a BC that uses no functions, e.g. critical depth conditions,
  // and are computed internally
  MultiPatch<double> crit_depth(createMultiPatchSpace(
    plist.sublist("critical depth"), mesh, Entity_kind::FACE, OPERATOR_BC_NEUMANN));

  for (int i = 0; i != crit_depth.size(); ++i) {
    auto ps = (*crit_depth.space)[i];
    auto ids = ps->getIDs();
    auto& patch_v = crit_depth[i].data;
    Kokkos::parallel_for(
      "crit_depth", ids.extent(0), KOKKOS_LAMBDA(const int& i) { patch_v(i, 0) = 3.14; });
  }

  // finally, one for a seepage-like condition, which may switch dynamically in
  // time between DIRICHLET and NEUMANN, and needs a function
  int OPERATOR_BC_CONDITIONAL = -1; // NOTE: -1 is necessary here to tell the
                                    // PatchSpace that it needs to store a VIEW of
                                    // flags, and not one flag
  MeshFunction mf_s(plist.sublist("seepage face"),
                    mesh,
                    "boundary pressure",
                    Entity_kind::FACE,
                    OPERATOR_BC_CONDITIONAL);
  MultiPatch<double> seepage(mf_s.createMPS(false));
  MultiPatch<int> seepage_flags(mf_s.createMPS(false));
  mf_s.Compute(1.0, seepage);
  const AmanziMesh::MeshCache& m = mesh->getCache();

  for (int i = 0; i != seepage.size(); ++i) {
    auto ps = (*seepage.space)[i];
    auto ids = ps->getIDs();
    auto& patch_v = seepage[i].data;
    auto& flags_v = seepage_flags[i].data;
    Kokkos::parallel_for(
      "seepage", ids.extent(0), KOKKOS_LAMBDA(const int& i) {
        auto fc = m.getFaceCentroid(ids(i));
        if (fc[2] > 2.0) {
          patch_v(i, 0) = 0.0;
          flags_v(i, 0) = OPERATOR_BC_NEUMANN;
        } else {
          flags_v(i, 0) = OPERATOR_BC_DIRICHLET;
        }
      });
  }

  // Operator::BCs objects consist of two CVs -- one for markers and one for values
  CompositeVectorSpace bc_space;
  bc_space.SetMesh(mesh)->AddComponent("face", Entity_kind::FACE, 1);
  auto cspace = bc_space.CreateSpace();
  CompositeVector_<int> bc_markers(cspace);
  CompositeVector bc_values(cspace);

  // copy in values
  copyMultiPatchToCompositeVector(dirichlet, "face", bc_values, bc_markers);
  copyMultiPatchToCompositeVector(neumann, "face", bc_values, bc_markers);
  copyMultiPatchToCompositeVector(crit_depth, "face", bc_values, bc_markers);
  copyMultiPatchToCompositeVector(seepage, "face", bc_values, bc_markers);

  // copy in flags
  copyMultiPatchToCompositeVector(seepage_flags, "face", bc_markers);

  CHECK_EQUAL(4, mesh->getSetSize("RIGHT", Entity_kind::FACE, Parallel_kind::OWNED));
  CHECK_EQUAL(4, mesh->getSetSize("LEFT", Entity_kind::FACE, Parallel_kind::OWNED));
  CHECK_EQUAL(4, mesh->getSetSize("TOP", Entity_kind::FACE, Parallel_kind::OWNED));
  CHECK_EQUAL(4, mesh->getSetSize("BOTTOM", Entity_kind::FACE, Parallel_kind::OWNED));
  CHECK_EQUAL(4, mesh->getSetSize("FRONT", Entity_kind::FACE, Parallel_kind::OWNED));
  CHECK_EQUAL(4, mesh->getSetSize("BACK", Entity_kind::FACE, Parallel_kind::OWNED));

  // check...
  {
    auto vals = bc_values.viewComponent<MemSpace_kind::HOST>("face", false);
    auto markers = bc_markers.viewComponent<MemSpace_kind::HOST>("face", false);

    for (int f = 0; f != vals.extent(0); ++f) {
      auto fc = mesh->getFaceCentroid(f);
      std::cout << "Checking " << f << " at " << fc << std::endl;
      if (fc[0] == 0.0) {
        // left
        CHECK_EQUAL(OPERATOR_BC_DIRICHLET, markers(f, 0));
        CHECK_EQUAL(1.0, vals(f, 0));
      } else if (fc[0] == 4.0) {
        // right
        CHECK_EQUAL(OPERATOR_BC_DIRICHLET, markers(f, 0));
        CHECK_EQUAL(1.0, vals(f, 0));
      } else if (fc[1] == 0.0) {
        // front
        CHECK_EQUAL(OPERATOR_BC_NEUMANN, markers(f, 0));
        CHECK_EQUAL(3.14, vals(f, 0));
      } else if (fc[1] == 4.0) {
        // back
        if (fc[2] > 2.0) {
          CHECK_EQUAL(OPERATOR_BC_NEUMANN, markers(f, 0));
          CHECK_EQUAL(0.0, vals(f, 0));
        } else {
          CHECK_EQUAL(OPERATOR_BC_DIRICHLET, markers(f, 0));
          CHECK_EQUAL(5.0, vals(f, 0));
        }
      } else if (fc[2] == 0.0) {
        // bottom
        CHECK_EQUAL(OPERATOR_BC_NEUMANN, markers(f, 0));
        CHECK_EQUAL(0.0, vals(f, 0));
      } else if (fc[2] == 4.0) {
        // top
        CHECK_EQUAL(OPERATOR_BC_NEUMANN, markers(f, 0));
        CHECK_EQUAL(3.0, vals(f, 0));
      } else {
        CHECK_EQUAL(0, markers(f, 0));
        CHECK_EQUAL(0.0, vals(f, 0));
      }
    }
  }
}
