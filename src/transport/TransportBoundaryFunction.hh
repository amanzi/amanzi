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
#include "MultiFunction.hh"
#include "unique_mesh_function.hh"

namespace Amanzi {
namespace AmanziTransport {

typedef std::pair<std::string, int> Action;

class TransportBoundaryFunction : public Functions::UniqueMeshFunction {
 public:
  TransportBoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh) :
      UniqueMeshFunction(mesh),
      finalized_(false) {};
  
  virtual void Compute(double time) = 0;

  // iterator methods
  typedef std::map<int,double>::const_iterator Iterator;
  Iterator begin() const { return value_.begin(); }
  Iterator end() const  { return value_.end(); }
  Iterator find(const int j) const { return value_.find(j); }
  std::map<int,double>::size_type size() { return value_.size(); }

 public:
  std::vector<std::string> tcc_name;

 protected:
  std::map<int,double> value_;
  bool finalized_;
};

}  // namespace AmanziTransport
}  // namespace Amanzi


#endif
