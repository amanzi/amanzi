/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      (v1) Daniil Svyatsky
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_DYNAMICBOUNDARY_FUNCTION_HH_
#define AMANZI_DYNAMICBOUNDARY_FUNCTION_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "MultiFunction.hh"
#include "UniqueMeshFunction.hh"
#include "BoundaryFunction.hh"

namespace Amanzi {
namespace Functions {

class DynamicBoundaryFunction : public BoundaryFunction {
 public:
  DynamicBoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : BoundaryFunction(mesh){};

  void AddFunction(const Teuchos::RCP<BoundaryFunction>& f);

  Teuchos::RCP<BoundaryFunction> GetFunction(int id) { return func_[id]; }

  int Func_ID(double time);

  void Compute(double time);

 protected:
  std::vector<Teuchos::RCP<BoundaryFunction>> func_;
};

} // namespace Functions
} // namespace Amanzi


#endif
