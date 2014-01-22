/* -------------------------------------------------------------------------
This is the Transport component of Amanzi

License: see $AMANZI_DIR/COPYRIGHT
Author (v1): Neil Carlson
       (v2): Ethan Coon

Function applied to a mesh component with at most one function 
application per entity.
------------------------------------------------------------------------- */

#ifndef AMANZI_TRANSPORT_BOUNDARY_FUNCTION_ALQUIMIA_HH_
#define AMANZI_TRANSPORT_BOUNDARY_FUNCTION_ALQUIMIA_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "MultiFunction.hh"
#include "TransportBoundaryFunction.hh"

namespace Amanzi {
namespace AmanziTransport {

class TransportBoundaryFunction_Alquimia : public TransportBoundaryFunction {
 public:
  TransportBoundaryFunction_Alquimia(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh) :
      TransportBoundaryFunction(mesh) {};
  
  void Compute(double time);

 private:
  void Define(const std::vector<std::string> &regions);

  void Define(std::string region);
};

}  // namespace AmanziTransport
}  // namespace Amanzi


#endif
