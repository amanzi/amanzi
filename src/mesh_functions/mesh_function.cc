#include "mesh_function.hh"

#include <algorithm>
#include "errors.hh"

namespace Amanzi {

/* ****************************************************************
* Similar to Define() but uses a single string.
**************************************************************** */
void MeshFunction::DefineFromString(const std::string region, const Teuchos::RCP<const Function>& f)
{
  std::vector<std::string> regions(1, region);
  Define(regions, f);
}

} // namespace Amanzi
