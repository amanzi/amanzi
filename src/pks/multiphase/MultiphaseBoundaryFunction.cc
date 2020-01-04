/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "MultiphaseBoundaryFunction.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor
****************************************************************** */
MultiphaseBoundaryFunction::MultiphaseBoundaryFunction(const Teuchos::ParameterList& plist)
{
  rainfall_ = false;
  if (plist.isParameter("rainfall")) 
    rainfall_ = plist.get<bool>("rainfall");
}


/* ****************************************************************
* Process additional parameters for BC submodels.
**************************************************************** */
void MultiphaseBoundaryFunction::ComputeSubmodel(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  int dim = mesh->space_dimension();

  if (rainfall_) {
    for (auto it = begin(); it != end(); ++it) {
      int f = it->first;
      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      it->second[0] *= fabs(normal[dim - 1]) / norm(normal);
    }
  }
}

}  // namespace Multiphase
}  // namespace Amanzi
