/*
  This is the Operator component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discrete diffusion operator of a surface.
*/

#ifndef AMANZI_OPERATOR_DIFFUSION_SURFACE_HH_
#define AMANZI_OPERATOR_DIFFUSION_SURFACE_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "tensor.hh"
#include "CompositeVector.hh"

#include "OperatorDiffusion.hh"
#include "OperatorTypeDefs.hh"
#include "NonlinearCoefficient.hh"

namespace Amanzi {
namespace Operators {

class OperatorDiffusionSurface : public OperatorDiffusion {
 public:
  OperatorDiffusionSurface() { InitDiffusionSurface_(); }
  OperatorDiffusionSurface(Teuchos::RCP<const CompositeVectorSpace> cvs, int dummy) 
      : OperatorDiffusion(cvs, dummy) { InitDiffusionSurface_(); }
  OperatorDiffusionSurface(const Operator& op) : OperatorDiffusion(op) { InitDiffusionSurface_(); };
  ~OperatorDiffusionSurface() {}; 

  // re-implementation of basic operator virtual members

 private:
  void InitDiffusionSurface_();
};

}  // namespace Operators
}  // namespace Amanzi


#endif


