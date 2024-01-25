/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  This smoothness indicator compares a reconstructed polynomial
  with that in neighbooring cells.
*/

#include <vector>

#include "Teuchos_RCP.hpp"

#include "MeshFramework.hh"

#include "Reconstruction.hh"
#include "SmoothnessIndicatorShu.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
*
****************************************************************** */
void
SmoothnessIndicatorShu::Init(Teuchos::ParameterList& plist)
{}


/* ******************************************************************
*
****************************************************************** */
void
SmoothnessIndicatorShu::Compute(const Teuchos::RCP<Reconstruction>& lifting)
{
  int ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  measure_.resize(ncells_owned);

  lifting->data()->scatterMasterToGhosted();

  AmanziMesh::Entity_ID_View cells;

  for (int c = 0; c < ncells_owned; ++c) {
    double value(1.0);
    const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
    auto poly = lifting->getPolynomial(c);

    cells = AmanziMesh::getCellFaceAdjacentCells(
      *mesh_, c, AmanziMesh::Parallel_kind::ALL);
    for (int c1 : cells) {
      auto poly1 = lifting->getPolynomial(c1);
      poly1.ChangeOrigin(xc);

      auto dpoly = poly1 - poly;
      double tmp = dpoly.normInf() / poly.normInf();
      value = std::min(value, 1.0 - tmp);
    }

    measure_[c] = value;
  }
}

} // namespace Operators
} // namespace Amanzi
