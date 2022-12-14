/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: (v1) Neil Carlson
           (v2) Ethan Coon
*/

/*
  Mesh Functions

  Function applied to a mesh component with at most one function
  application per entity.
*/

#include "DynamicBoundaryFunction.hh"

namespace Amanzi {
namespace Functions {

void
DynamicBoundaryFunction::AddFunction(const Teuchos::RCP<BoundaryFunction>& f)
{
  func_.push_back(f);
}

int
DynamicBoundaryFunction::Func_ID(double time)
{
  // lazily generate space for the values
  if (!finalized_) { Finalize(); }

  if (unique_specs_.size() == 0) return 0;

  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);
  args[0] = time;

  UniqueSpecList::const_iterator uspec = unique_specs_[AmanziMesh::FACE]->begin();
  int val = std::floor((*(*uspec)->first->second)(args)[0]);

  return val;
}


// Evaluate values at time.
void
DynamicBoundaryFunction::Compute(double time)
{
  // lazily generate space for the values
  if (!finalized_) { Finalize(); }

  if (unique_specs_.size() == 0) return;

  int func_id = Func_ID(time);

  if (func_id >= func_.size()) return;
  func_[func_id]->Compute(time);
};

} // namespace Functions
} // namespace Amanzi
