/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#include <cmath>
#include <string>
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "Mesh.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Flow;

/* ******************************************************************
* Calculate L2 error in pressure.
****************************************************************** */
double
CalculatePressureCellError(Teuchos::RCP<const Mesh> mesh, const Epetra_MultiVector& p)
{
  double k1 = 0.5, k2 = 2.0, g = 2.0, a = 5.0, cr = 1.02160895462971866; // analytical data
  double f1 = sqrt(1.0 - g * k1 / cr);
  double f2 = sqrt(g * k2 / cr - 1.0);

  int dim = mesh->getSpaceDimension();
  double pexact, error_L2 = 0.0;
  for (int c = 0; c < p.MyLength(); c++) {
    const AmanziGeometry::Point& xc = mesh->getCellCentroid(c);
    double volume = mesh->getCellVolume(c);

    double z = xc[dim - 1];
    if (z < -a) {
      pexact = f1 * tan(cr * (z + 2 * a) * f1 / k1);
    } else {
      pexact = -f2 * tanh(cr * f2 * (z + a) / k2 - atanh(f1 / f2 * tan(cr * a * f1 / k1)));
    }
    error_L2 += std::pow(p[0][c] - pexact, 2.0) * volume;
    // std::cout << z << " " << p[0][c] << " exact=" <<  pexact << std::endl;
  }
  return sqrt(error_L2);
}


/* ******************************************************************
* Calculate l2 error (small l) in darcy flux.
****************************************************************** */
double
CalculateDarcyFluxError(Teuchos::RCP<const Mesh> mesh, const Epetra_MultiVector& flux)
{
  int dim = mesh->getSpaceDimension();
  AmanziGeometry::Point velocity_exact(dim);

  double cr = 1.02160895462971866; // analytical data
  velocity_exact[dim - 1] = -cr;

  int nfaces_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);

  double error_l2 = 0.0;
  for (int f = 0; f < nfaces_owned; f++) {
    const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);
    error_l2 += std::pow(flux[0][f] - velocity_exact * normal, 2.0);
    // std::cout << f << " " << flux[0][f] << " exact=" << velocity_exact * normal
    //           << " xf=" << mesh->getFaceCentroid(f) << std::endl;
  }
  return sqrt(error_l2 / nfaces_owned);
}


/* ******************************************************************
* Calculate L2 divergence error in darcy flux.
****************************************************************** */
double
CalculateDarcyDivergenceError(Teuchos::RCP<const Mesh> mesh, const Epetra_MultiVector& flux)
{
  double error_L2 = 0.0;
  int ncells_owned =
    mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;

    mesh->getCellFacesAndDirections(c, &faces, &dirs);
    int nfaces = faces.size();

    double div = 0.0;
    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      div += flux[0][f] * dirs[i];
    }
    error_L2 += div * div / mesh->getCellVolume(c);
    //std::cout << c << " div=" << div << " exact=0.0" << std::endl;
  }
  return sqrt(error_L2);
}
