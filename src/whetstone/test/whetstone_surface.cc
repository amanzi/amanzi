/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "Point.hh"
#include "SingleFaceMesh.hh"

#include "MFD3D_Diffusion.hh"
#include "Tensor.hh"


/* **************************************************************** */
TEST(DARCY_SURFACE_MESH)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Mass matrix for surface mesh" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 3, 4);

 Tensor<> T(2, 1);
  T(0, 0) = 1;

  for (int f = 0; f < mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL); ++f) {
    RCP<AmanziMesh::SingleFaceMesh> surfmesh = Teuchos::rcp(new AmanziMesh::SingleFaceMesh(mesh, f));

    MFD3D_Diffusion mfd(surfmesh);

    DenseMatrix<> M;
    mfd.MassMatrix(0, T, M);
    int nrows = M.NumRows();
    
    for (int i = 0; i < nrows; ++i) {
      for (int j = 0; j < nrows; ++j) {
        double val = (i != j) ? 0.0 : surfmesh->getCellVolume(0) / 2;
        CHECK_CLOSE(M(i, j), val, 1e-10);
      }
    }
  }
}


/* **************************************************************** */
TEST(DARCY_SURFACE)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "\nTest: Mass matrix for cell face" << std::endl;
#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  auto comm = Amanzi::getCommSelf();
#endif

  MeshFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  RCP<Mesh> mesh = meshfactory.create("test/surface.exo");

  MFD3D_Diffusion mfd(mesh);

 Tensor<> T(2, 1);
  T(0, 0) = 1;

  for (int c = 0; c < 3; c++) {
    auto faces = mesh->getCellFaces(c);
    int nfaces = faces.size();

    DenseMatrix<> W(nfaces, nfaces);
    for (int method = 0; method < 2; method++) {
      int ok(0);
      if (method == 0) {
        ok = mfd.MassMatrixInverseSurface(c, T, W);
      } else if (method == 1) {
        if (c == 2) continue;
        ok = mfd.MassMatrixInverseSurfaceMMatrix(c, T, W);
        std::cout << "Number of simplex itrs=" << mfd.simplex_num_itrs() << " code=" << ok
                  << std::endl;
        std::cout << "Functional value=" << mfd.simplex_functional() << std::endl;
      }

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
      double vyy = 0.0, vxy = 0.0, volume = mesh->getCellVolume(c);
      for (int i = 0; i < nfaces; i++) {
        int f = faces[i];
        yi = mesh->getFaceNormal(f, c)[1]; 
        for (int j = 0; j < nfaces; j++) {
          f = faces[j];
          xj = mesh->getFaceNormal(f, c)[0];
          yj = mesh->getFaceNormal(f, c)[1];
          vxy += W(i, j) * yi * xj;
          vyy += W(i, j) * yi * yj;
        }
      }
      CHECK_CLOSE(vyy, volume, 1e-10);
      CHECK_CLOSE(vxy, 0.0, 1e-10);
    }
  }
}
