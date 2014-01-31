/* -------------------------------------------------------------------------
This is the Transport component of Amanzi

License: see $AMANZI_DIR/COPYRIGHT
Author (v1): Neil Carlson
       (v2): Ethan Coon
       (v3): Konstantin Lipnikov

Function applied to a mesh component with at most one function 
application per entity.
------------------------------------------------------------------------- */

#ifndef AMANZI_TRANSPORT_BOUNDARY_FUNCTION_HH_
#define AMANZI_TRANSPORT_BOUNDARY_FUNCTION_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "unique_mesh_function.hh"

namespace Amanzi {
namespace AmanziTransport {

typedef std::pair<std::string, int> Action;

class TransportBoundaryFunction : public Functions::UniqueMeshFunction {
 public:
  TransportBoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh) :
      UniqueMeshFunction(mesh) {};
  
  virtual void Compute(double time) = 0;

  // access
  std::vector<string>& tcc_names() { return tcc_names_; }
  std::vector<int>& tcc_index() { return tcc_index_; }

  std::vector<int>& faces() { return faces_; }
  std::vector<std::vector<double> >& values() { return values_; }

 protected:
  std::vector<std::string> tcc_names_;  // list of component names
  std::vector<int> tcc_index_;  // index of component in the global list

  std::vector<int> faces_;  // list of boundary faces 
  std::vector<std::vector<double> > values_;  // component values on boundary faces
};

}  // namespace AmanziTransport
}  // namespace Amanzi


#endif
