/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Navier Stokes PK

*/

#include "WhetStoneDefs.hh"

#include "NavierStokesBoundaryFunction.hh"

namespace Amanzi {
namespace NavierStokes {

/* ****************************************************************
* Constructor: extract attributes to setup a submodel.
**************************************************************** */
NavierStokesBoundaryFunction::NavierStokesBoundaryFunction(const Teuchos::ParameterList& plist) {}


/* ****************************************************************
* Process additional parameters for BC submodels.
**************************************************************** */
void
NavierStokesBoundaryFunction::ComputeSubmodel(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  int dim = mesh->getSpaceDimension();

  if (type_ == WhetStone::DOF_Type::NORMAL_COMPONENT) {
    for (auto it = begin(); it != end(); ++it) {
      int f = it->first;
      const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);

      double tmp(0.0);
      for (int k = 0; k < dim; ++k) { tmp += it->second[k] * normal[k]; }
      tmp /= norm(normal);
      it->second = std::vector<double>(1, tmp);
    }
  }
}

} // namespace NavierStokes
} // namespace Amanzi
