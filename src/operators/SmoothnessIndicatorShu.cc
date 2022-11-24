/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  This smoothness indicator compares a reconstructed polynomial
  with that in neighbooring cells.
*/

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

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
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  measure_.resize(ncells_owned);

  lifting->data()->ScatterMasterToGhosted();

  AmanziMesh::Entity_ID_List cells;

  for (int c = 0; c < ncells_owned; ++c) {
    double value(1.0);
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    auto poly = lifting->getPolynomial(c);

    mesh_->cell_get_face_adj_cells(c, AmanziMesh::Parallel_type::ALL, &cells);
    for (int c1 : cells) {
      auto poly1 = lifting->getPolynomial(c1);
      poly1.ChangeOrigin(xc);

      auto dpoly = poly1 - poly;
      double tmp = dpoly.NormInf() / poly.NormInf();
      value = std::min(value, 1.0 - tmp);
    }

    measure_[c] = value;
  }
}

} // namespace Operators
} // namespace Amanzi
