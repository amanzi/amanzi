#ifndef AMANZI_DOMAIN_FUNCTION_HH_
#define AMANZI_DOMAIN_FUNCTION_HH_

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "function.hh"
#include "mesh_function.hh"

namespace Amanzi {

class DomainFunction : public MeshFunction {
 public:
  DomainFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) { mesh_ = mesh; }

  void Define(const std::vector<std::string> &regions, const Teuchos::RCP<const Function>& f);
};

} // namespace Amanzi

#endif  // AMANZI_DOMAIN_FUNCTION_HH_
