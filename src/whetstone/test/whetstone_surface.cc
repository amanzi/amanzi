/*
  The discretization component of Amanzi.
  License: BSD
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "MeshFactory.hh"
#include "MeshAudit.hh"

#include "AmanziComm.hh"
#include "Mesh.hh"
#include "Point.hh"

#include "MFD3D_Diffusion.hh"
#include "Tensor.hh"


/* **************************************************************** */
TEST(DARCY_SURFACE)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Mass matrix for surface Darcy" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh = meshfactory.create("test/surface.exo");

  MFD3D_Diffusion mfd(mesh);

  Tensor T(2, 1);
  T(0, 0) = 1;

  for (int c = 0; c < 3; c++) {
    Kokkos::View<Amanzi::WhetStone::Entity_ID*> faces;
    Kokkos::View<int*> dirs;

    mesh->cell_get_faces_and_dirs(c, faces, dirs);
    int nfaces = faces.extent(0);

    DenseMatrix W(nfaces, nfaces);
    int ok = mfd.MassMatrixInverseSurface(c, T, W);

    printf("Inverse mass matrix for cell %d  err=%d\n", c, ok);
    for (int i = 0; i < nfaces; i++) {
      for (int j = 0; j < nfaces; j++) printf("%8.4f ", W(i, j));
      printf("\n");
    }

    // verify SPD propery
    for (int i = 0; i < nfaces; i++) CHECK(W(i, i) > 0.0);

    // verify exact integration property
    W.Inverse();

    double xj, yi, yj;
    double vyy = 0.0, vxy = 0.0, volume = mesh->cell_volume(c, false);
    for (int i = 0; i < nfaces; i++) {
      int f = faces(i);
      yi = mesh->face_normal(f)[1] * dirs(i);
      for (int j = 0; j < nfaces; j++) {
        f = faces(j);
        xj = mesh->face_normal(f)[0] * dirs(j);
        yj = mesh->face_normal(f)[1] * dirs(j);
        vxy += W(i, j) * yi * xj;
        vyy += W(i, j) * yi * yj;
      }
    }
    CHECK_CLOSE(vyy, volume, 1e-10);
    CHECK_CLOSE(vxy, 0.0, 1e-10);
  }
}
