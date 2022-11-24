/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <iterator>

#include "MultiphaseBoundaryFunction.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor defaults to liquid component which is water
****************************************************************** */
MultiphaseBoundaryFunction::MultiphaseBoundaryFunction(const Teuchos::ParameterList& plist)
{
  rainfall_ = false;
  if (plist.isParameter("name")) component_name_ = plist.get<std::string>("name");

  component_phase_ = MULTIPHASE_PHASE_LIQUID;
  if (plist.isParameter("phase")) {
    std::string phase = plist.get<std::string>("phase");
    if (phase == "gas") component_phase_ = MULTIPHASE_PHASE_GAS;
    if (phase == "napl") component_phase_ = MULTIPHASE_PHASE_NAPL;
  }

  if (plist.isParameter("rainfall")) rainfall_ = plist.get<bool>("rainfall");
}


/* ****************************************************************
* Process additional parameters for BC submodels.
**************************************************************** */
void
MultiphaseBoundaryFunction::ComputeSubmodel(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
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


/* ****************************************************************
* Find position of component in the list of names
**************************************************************** */
void
MultiphaseBoundaryFunction::SetComponentId(const std::vector<std::string>& names)
{
  auto it = std::find(names.begin(), names.end(), component_name_);
  component_id_ = (it == names.end()) ? -1 : std::distance(names.begin(), it);
}

} // namespace Multiphase
} // namespace Amanzi
