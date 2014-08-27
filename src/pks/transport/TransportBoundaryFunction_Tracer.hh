/* -------------------------------------------------------------------------
This is the Transport component of Amanzi

License: see $AMANZI_DIR/COPYRIGHT
Author (v1): Neil Carlson
       (v2): Ethan Coon
       (v3): Konstantin Lipnikov

Function applied to a mesh component with at most one function 
application per entity.
------------------------------------------------------------------------- */

#ifndef AMANZI_TRANSPORT_BOUNDARY_FUNCTION_TRACER_HH_
#define AMANZI_TRANSPORT_BOUNDARY_FUNCTION_TRACER_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "MultiFunction.hh"
#include "TransportBoundaryFunction.hh"

namespace Amanzi {
namespace Transport {

class TransportBoundaryFunction_Tracer : public TransportBoundaryFunction {
 public:
  TransportBoundaryFunction_Tracer(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh) :
      TransportBoundaryFunction(mesh) {};
  ~TransportBoundaryFunction_Tracer() {};
  
  void Compute(double time);

  void Define(const std::vector<std::string> &regions,
              const Teuchos::RCP<const MultiFunction> &f);

  void Define(std::string region,
              const Teuchos::RCP<const MultiFunction> &f);

 private:
  void Finalize_();
};

}  // namespace Transport
}  // namespace Amanzi


#endif
