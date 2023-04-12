/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>
#include <vector>

#include "UnitTest++.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"

#include "MeshFactory.hh"
#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "AmanziVector.hh"
#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;


struct TestHarness {
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;

  TestHarness()
  {
    comm = getDefaultComm();
    MeshFactory meshfactory(comm);
    mesh = meshfactory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 8, 1, 1);
  }

  ~TestHarness() {}

  Teuchos::RCP<CompositeVector>
  CreateVec(const std::string& name, Entity_kind kind, int num_vectors)
  {
    auto x_space = Teuchos::rcp(new CompositeVectorSpace());
    x_space->SetMesh(mesh)->SetGhosted()->SetComponents(
      {
        name,
      },
      {
        kind,
      },
      {
        num_vectors,
      });
    return x_space->Create();
  }
};


SUITE(COMMON_MESH_OPERATIONS)
{
  TEST_FIXTURE(TestHarness, FOR_EACH_CELL_VOLUME_LAMBDA)
  {
    // does a Kokkos for_each over cell volumes to calculate a derived quantity
    // using lambdas
    auto sl = CreateVec("cell", Entity_kind::CELL, 1);
    sl->putScalar(0.5);
    Teuchos::RCP<const CompositeVector> sl_c(sl);

    auto poro = CreateVec("cell", Entity_kind::CELL, 1);
    poro->putScalar(0.25);
    Teuchos::RCP<const CompositeVector> poro_c(poro);

    auto wc = CreateVec("cell", Entity_kind::CELL, 1);

    // compute on device
    {
      auto sl_view = sl_c->viewComponent<DefaultDevice>("cell", 0, false);
      auto poro_view = poro_c->viewComponent<DefaultDevice>("cell", 0, false);
      auto wc_view = wc->viewComponent<DefaultDevice>("cell", 0, false);
      Mesh* m = mesh.get();

      // typedef Kokkos::RangePolicy<DefaultDevice, LO> Policy_type;
      Kokkos::parallel_for(
        "data_structures_common_patterns::FOR_EACH_CELL_VOLUME_LAMBDA",
        sl_view.extent(0),
        KOKKOS_LAMBDA(const LO i) {
          wc_view[i] = sl_view[i] * poro_view[i] * m->getCellVolume(i);
        });
    }

    // test on the host
    {
      Teuchos::RCP<const CompositeVector> wc_c(wc);
      auto wc_view = wc_c->viewComponent<Amanzi::HostSpaceSpecial>("cell", 0, false);
      for (LO i = 0; i != wc_view.extent(0); ++i) {
        CHECK_CLOSE(0.5 * 0.25 * 8.0, wc_view(i), 1.e-10);
      }
    }
  }


  TEST_FIXTURE(TestHarness, CALC_TRANSMISSIBILITY)
  {
    // Computes transmissibility through flat parallelism
    auto trans = CreateVec("face", Entity_kind::FACE, 1);
    {
      auto trans_view = trans->viewComponent<DefaultDevice>("face", 0, false);
      Mesh* m = mesh.get();

      typedef Kokkos::RangePolicy<DefaultDevice, LO> Policy_type;
      Kokkos::parallel_for(
        "data_structures_common_patterns::CALC_TRANSMISSIBILITY",
        trans_view.extent(0),
        KOKKOS_LAMBDA(const LO& f) {
          // face info
          const auto& fc = m->getFaceCentroid(f);
          const auto& f_normal = m->getFaceNormal(f);
          double area = m->getFaceArea(f);

          // neighbor cells
          auto cellids = m->getFaceCells(f);

          // this could be a kokkos reduction...
          for (int i = 0; i != cellids.size(); ++i) {
            const auto& cc = m->getCellCentroid(cellids[i]);
            auto bisector = fc - cc;
            double s = area / norm(bisector);

            // cancellation here, but not in the real
            // usage, so please don't remove the
            // cancellation...
            double perm = (bisector * f_normal) * s;
            double dxn = bisector * f_normal;
            trans_view(f) += fabs(dxn / perm);
          }
          trans_view(f) = 1.0 / trans_view(f);
        });
    }
  }
}
