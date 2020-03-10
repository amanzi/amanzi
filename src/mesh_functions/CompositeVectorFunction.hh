/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
    Ethan Coon (coonet@ornl.gov)
*/

//!

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "MeshFunction.hh"
#include "CompositeVector.hh"

namespace Amanzi {
namespace Functions {

class CompositeVectorFunction {
 public:
  CompositeVectorFunction(const MultiPatchSpace& space,
                          const std::vector<Teuchos::RCP<const MultiFunction>>& funcs) :
      space_(space),
      funcs_(funcs) {}    

  void set_mesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) {
    space_.set_mesh(mesh);
    
  }
  
  void Compute(double time, CompositeVector& vec) {
    Functions::computeMeshFunction(funcs_, time, space_, vec);
  }

 protected:
  MultiPatchSpace space_;
  std::vector<Teuchos::RCP<const MultiFunction>> funcs_;
};


//
// Creates a function without a mesh, which must be set later.  For use by
// evaluators.
//
inline Teuchos::RCP<CompositeVectorFunction>
createCompositeVectorFunction(Teuchos::ParameterList& plist) {
  auto speclist_funcs = processListWithFunction(plist);
  return Teuchos::rcp(new CompositeVectorFunction(speclist_funcs.first,
          speclist_funcs.second));
}

//
// Creates a function with a mesh.  Preferred version of the above, which
// should get deprecated eventually.
//
inline Teuchos::RCP<CompositeVectorFunction>
createCompositeVectorFunction(Teuchos::ParameterList& plist,
        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) {
  auto func = createCompositeVectorFunction(plist);
  func->set_mesh(mesh);
  return func;
}



} // namespace Functions
} // namespace Amanzi


