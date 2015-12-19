/*
  Energy PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Function applied to a mesh component with at most one function 
  application per entity. So far, this class delegates work to 
  the base class.
*/

#ifndef AMANZI_ENERGY_BOUNDARY_FUNCTION_HH_
#define AMANZI_ENERGY_BOUNDARY_FUNCTION_HH_

#include <map>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CommonDefs.hh"
#include "Mesh.hh"
#include "MultiFunction.hh"
#include "PK_BoundaryFunction.hh"
#include "UniqueMeshFunction.hh"

namespace Amanzi {
namespace Energy {

typedef std::pair<std::string, int> Action;

class EnergyBoundaryFunction : public PK_BoundaryFunction {
 public:
  EnergyBoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PK_BoundaryFunction(mesh) {};
  ~EnergyBoundaryFunction() {};
};

}  // namespace Energy
}  // namespace Amanzi

#endif
