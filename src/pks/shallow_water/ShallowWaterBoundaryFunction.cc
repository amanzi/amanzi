/*
 Shallow water PK
 
 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.
 
 Author: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "WhetStoneDefs.hh"

#include "ShallowWaterBoundaryFunction.hh"

namespace Amanzi {
namespace ShallowWater {
        
/* ****************************************************************
* Constructor: extract attributes to setup a submodel.
**************************************************************** */
ShallowWaterBoundaryFunction::ShallowWaterBoundaryFunction(const Teuchos::ParameterList& plist)
{
}


/* ****************************************************************
* Process additional parameters for BC submodels.
**************************************************************** */
void ShallowWaterBoundaryFunction::ComputeSubmodel(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  int dim = mesh->space_dimension();

  if (type_ == WhetStone::DOF_Type::NORMAL_COMPONENT) {
    for (auto it = begin(); it != end(); ++it) {
      int f = it->first;
      const AmanziGeometry::Point& normal = mesh->face_normal(f);

      double tmp(0.0);
      for (int k = 0; k < dim; ++k) {
        tmp += it->second[k] * normal[k];
      }
      tmp /= norm(normal);
      it->second = std::vector<double>(1, tmp);
    }
  }
}
        
}  // namespace ShallowWater
}  // namespace Amanzi

