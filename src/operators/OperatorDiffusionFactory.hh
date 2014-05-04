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

class OperatorDiffusionFactory {
 public:
  OperatorDiffusionFactory() {};
  ~OperatorDiffusionFactory() {};

  Teuchos::RCP<OperatorDiffusion> Create(Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                                         const Teuchos::ParameterList& op_list);

  Teuchos::RCP<OperatorDiffusion> Create(Teuchos::RCP<Operator> op,
                                         const Teuchos::ParameterList& op_list);
};

}  // namespace Operators
}  // namespace Amanzi

#endif
