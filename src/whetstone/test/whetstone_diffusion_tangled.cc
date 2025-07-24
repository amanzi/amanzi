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
#include <utility>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

#include "MeshFrameworkFactory.hh"
#include "Point.hh"

#include "MFD3D_Diffusion.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace AmanziMesh {

struct MeshAlgorithmsTangled : public MeshAlgorithms {
  virtual ~MeshAlgorithmsTangled() = default;

  virtual std::pair<double, AmanziGeometry::Point> computeCellGeometry(const MeshHost& mesh,
                                                                       const Entity_ID c) const
  {
    double area = 0.0;
    AmanziGeometry::Point centroid(2);

    auto coords = mesh.getCellCoordinates(c);
    std::size_t np = coords.size();

    for (int i = 0; i < np; i++) centroid += coords[i];
    centroid /= np;

    // Compute the area as the sum of triangle areas
    for (int i = 0; i < np; i++) {
      AmanziGeometry::Point v1 = coords[i] - centroid;
      AmanziGeometry::Point v2 = coords[(i + 1) % np] - centroid;
      AmanziGeometry::Point v3 = 0.5 * v1 ^ v2;
      area += v3[0];
    }

    return std::make_pair(area, centroid);
  }
};

} // namespace AmanziMesh
} // namespace Amanzi


/* **************************************************************** */
TEST(INVERTED_CELL_2D)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: Mass matrix for inverted cell in 2D" << std::endl;
  auto comm = Amanzi::getDefaultComm();

  MeshFrameworkFactory meshfactory(comm);
  meshfactory.set_preference(Preference({ Framework::MSTK }));
  auto mesh_fw = meshfactory.create(0.0, 0.0, 1.1, 1.0, 1, 1);

  auto mesh = Teuchos::rcp(new AmanziMesh::MeshCache<MemSpace_kind::HOST>(
    mesh_fw, Teuchos::rcp(new MeshAlgorithmsTangled()), Teuchos::null));

  AmanziGeometry::Point x3(-1.5, 1.0);
  mesh->setNodeCoordinate(3, x3);
  mesh->recacheGeometry();
  std::cout << "oriented volume:   " << mesh->getCellVolume(0) << std::endl;
  std::cout << "modified centroid: " << mesh->getCellCentroid(0) << std::endl;

  MFD3D_Diffusion mfd(mesh);
  DeRham_Face drc(mfd);

  int nfaces = 4, cell = 0;
  DenseMatrix M;

  Tensor T(2, 2);
  T(0, 0) = 2.0;
  T(1, 1) = 2.0;
  T(0, 1) = 0.5;
  T(1, 0) = 0.5;

  M.Reshape(nfaces, nfaces);
  mfd.MassMatrix(cell, T, M);

  printf("Mass matrix for cell %3d\n", cell);
  PrintMatrix(M, "%8.4f ");

  // verify exact integration property
  auto [faces, dirs] = mesh->getCellFacesAndDirections(cell);

  double xi, yi, xj;
  double vxx = 0.0, vxy = 0.0, volume = mesh->getCellVolume(cell);
  for (int i = 0; i < nfaces; i++) {
    int f1 = faces[i];
    for (int j = 0; j < nfaces; j++) {
      int f2 = faces[j];

      xi = mesh->getFaceNormal(f1)[0] * dirs[i];
      yi = mesh->getFaceNormal(f1)[1] * dirs[i];
      xj = mesh->getFaceNormal(f2)[0] * dirs[j];

      xi /= mesh->getFaceArea(f1);
      yi /= mesh->getFaceArea(f1);
      xj /= mesh->getFaceArea(f2);

      vxx += M(i, j) * xi * xj;
      vxy += M(i, j) * yi * xj;
    }
  }

  CHECK_CLOSE(T(0, 0) * volume, vxx, 1e-10);
  CHECK_CLOSE(T(1, 0) * volume, vxy, 1e-10);
}
