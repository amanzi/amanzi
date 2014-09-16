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

namespace Amanzi {
namespace Operators {

class BCs;

class OperatorDiffusionSurface : public OperatorDiffusion {
 public:
  OperatorDiffusionSurface() {};
  OperatorDiffusionSurface(Teuchos::RCP<const CompositeVectorSpace> cvs, 
                           Teuchos::ParameterList& plist, Teuchos::RCP<BCs> bc) 
      : OperatorDiffusion(cvs, plist, bc) { InitDiffusionSurface_(plist); }
  OperatorDiffusionSurface(const Operator& op, 
                           Teuchos::ParameterList& plist, Teuchos::RCP<BCs> bc)
      : OperatorDiffusion(op, plist, bc) { InitDiffusionSurface_(plist); };
  ~OperatorDiffusionSurface() {}; 

  // re-implementation of basic operator virtual members

 private:
  void InitDiffusionSurface_(Teuchos::ParameterList& plist);
};

}  // namespace Operators
}  // namespace Amanzi


#endif


