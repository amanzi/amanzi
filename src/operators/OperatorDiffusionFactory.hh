/*
  This is the operators component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_FACTORY_HH_
#define AMANZI_OPERATOR_DIFFUSION_FACTORY_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "OperatorDiffusion.hh"

namespace Amanzi {
namespace Operators {

class BCs;

class OperatorDiffusionFactory {
 public:
  OperatorDiffusionFactory() {};
  ~OperatorDiffusionFactory() {};

  Teuchos::RCP<OperatorDiffusion> Create(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                                         Teuchos::RCP<BCs> bc,
                                         Teuchos::ParameterList& oplist,
                                         const AmanziGeometry::Point& g,
                                         int upwind_method);

  Teuchos::RCP<OperatorDiffusion> Create(Teuchos::RCP<Operator> op,
                                         Teuchos::RCP<BCs> bc,
                                         Teuchos::ParameterList& oplist);
};

}  // namespace Operators
}  // namespace Amanzi

#endif
