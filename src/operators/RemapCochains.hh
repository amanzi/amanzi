/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Remap of 2-cochains and 3-cochains with supporting tools
  for physics packages.

  The meshes are assumed to be topolgically identical.
*/

#ifndef AMANZI_OPERATOR_REMAP_COCHAINS_HH_
#define AMANZI_OPERATOR_REMAP_COCHAINS_HH_

#include <utility>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "Point.hh"
#include "CompositeVector.hh"
#include "dbc.hh"
#include "Point.hh"

namespace Amanzi {
namespace Operators {

class RemapCochains {
 public:
  RemapCochains(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
                const Teuchos::RCP<const AmanziMesh::Mesh> mesh1);
  ~RemapCochains() {};

  template<int K>
  void LowOrder(const CompositeVector& in, CompositeVector& out);
  template<int K>
  void HighOrder(const CompositeVector& in, CompositeVector& out);

  AmanziGeometry::Point getSweptEdgeCentroid(int e) const;
  std::vector<std::pair<AmanziGeometry::Point, AmanziGeometry::Point>>
  getSweptEdgeQuadrature(int e) const;
  AmanziGeometry::Point getSweptFaceCentroid(int f) const { return AmanziGeometry::Point(3); }

  // reconstruction
  Teuchos::RCP<CompositeVector>
  getCellFieldFromEdges(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, const CompositeVector& in);
  Teuchos::RCP<CompositeVector>
  getCellFieldFromFaces(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, const CompositeVector& in);

 private:
  int d_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh0_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh1_;
};

} // namespace Operators
} // namespace Amanzi

#endif
