#ifndef AMANZI_BOUNDARY_FUNCTION_HH_
#define AMANZI_BOUNDARY_FUNCTION_HH_

#include <vector>
#include <string>
#include <utility>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "function.hh"
#include "mesh_function.hh"

namespace Amanzi {

const int BOUNDARY_FUNCTION_ACTION_NONE = 0;
const int BOUNDARY_FUNCTION_ACTION_HEAD_RELATIVE = 1;

typedef std::pair<std::string, int> Action;

class BoundaryFunction : public MeshFunction {
 public:
  BoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) { mesh_ = mesh; }

  void Define(const std::vector<std::string>& regions,
              const Teuchos::RCP<const Function>& f, int action);
  void Compute(double T);
  void ComputeShift(double T, double* shift);

  const std::vector<Action>& actions() { return actions_; } 

 private:
  std::vector<Action> actions_;
};

} // namespace Amanzi

#endif  // AMANZI_BOUNDARY_FUNCTION_HH_
